package org.broadinstitute.hellbender.utils.minimap2;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/** Controls a minimap2 index file.
 * Creating an instance causes an entire index file (which may be many gigabytes) to be read.
 * It's an expensive operation that should be minimized.
 * It should be closed when you're done with it to reclaim large amounts of (non-Java) memory.
 *
 * This class is thread-safe, and should be shared among threads to conserve resources.
 *
 * Current implementation only reads the first chunk of multi-chunk indices.
 * (However, an entire human reference genome is typically in one chunk, so this is probably not a big issue.)
 */
public class MiniMap2Index implements AutoCloseable {
    private long nativeAddress;
    private int refCount;
    private List<String> refNames = null;

    private static boolean nativeLibLoaded = false;

    public MiniMap2Index( final String mmiFile ) {
        loadNativeLibrary();
        nativeAddress = openIndex(mmiFile);
        if ( nativeAddress == 0L ) {
            throw new MiniMap2Exception("Can't open minimap2 index file: " + mmiFile);
        }
        refCount = 0;
    }

    public boolean isOpen() { return nativeAddress != 0; }

    public void close() {
        final long addr;
        synchronized (this) {
            if ( (addr = nativeAddress) == 0L ) return;
            if ( refCount > 0 ) {
                throw new MiniMap2Exception("Can't close index:  it's in use.");
            }
            nativeAddress = 0L;
        }
        destroyIndex(addr);
    }

    /** returns an immutable list of contig names for the reference */
    public List<String> getRefNames() {
        if ( refNames == null ) {
            try {
                final long addr;
                synchronized (this) {
                    refCount += 1;
                    addr = nativeAddress;
                }
                if ( addr == 0L ) {
                    throw new MiniMap2Exception("Can't get ref names:  index is closed.");
                }
                final ByteBuffer refNameBuffer = getRefNames(addr);
                if ( refNameBuffer == null ) {
                    throw new MiniMap2Exception("Couldn't retrieve ref names:  Reasons are obscure.");
                }
                final byte[] nameBytes = new byte[refNameBuffer.capacity()];
                refNameBuffer.get(nameBytes, 0, nameBytes.length);
                destroyByteBuffer(refNameBuffer);
                final List<String> names = new ArrayList<>();
                int startIdx = 0;
                for ( int idx = 0; idx != nameBytes.length; ++idx ) {
                    if ( nameBytes[idx] == 0 ) {
                        names.add(new String(nameBytes, startIdx, idx - startIdx));
                        startIdx = idx + 1;
                    }
                }
                refNames = Collections.unmodifiableList(names);
            } finally {
                synchronized (this) {
                    refCount -= 1;
                }
            }
        }
        return refNames;
    }

    /** get the contig name for a MiniMap2Alignment's refId */
    public String getRefName( final int refId ) {
        return getRefNames().get(refId);
    }

    ByteBuffer createOptions( final String preset ) {
        synchronized (this) {
            if ( nativeAddress == 0L ) {
                throw new MiniMap2Exception("Can't create options for aligner:  index is closed.");
            }
            return createOptions(nativeAddress, preset);
        }
    }

    ByteBuffer align( final ByteBuffer opts, final ByteBuffer seqs ) {
        try {
            final long addr;
            synchronized (this) {
                refCount += 1;
                addr = nativeAddress;
            }
            if ( addr == 0L ) {
                throw new MiniMap2Exception("Can't align:  index is closed.");
            }
            return createAlignments(addr, opts, seqs);
        } finally {
            synchronized (this) {
                refCount -= 1;
            }
        }
    }

    private static synchronized void loadNativeLibrary() {
        if ( nativeLibLoaded ) return;

        final String libNameOverride = System.getProperty("LIBMM2_PATH");
        if ( libNameOverride != null ) {
            System.load(libNameOverride);
            nativeLibLoaded = true;
            return;
        }

        final String osName = System.getProperty("os.name", "unknown").toUpperCase();
        final String osArch = System.getProperty("os.arch");
        if ( !"x86_64".equals(osArch) && !"amd64".equals(osArch) ) {
            throw new MiniMap2Exception(
                    "We have pre-built minimap2 native libraries only for x86_64 and amd64 architectures. "+
                            "Your os.arch is " + osArch + ". " +
                            "Set property LIBMM2_PATH to point to a native library for your architecture.");
        }

        final String libName;
        if ( osName.startsWith("MAC") ) libName = "/libmm2.Darwin.dylib";
        else if ( osName.startsWith("LINUX") ) libName = "/libmm2.Linux.so";
        else {
            throw new MiniMap2Exception(
                    "We have pre-built minimap2 native libraries only for Linux and Mac. "+
                            "Your os.name is " + osName + ". " +
                            "Set property LIBMM2_PATH to point to a native library for your operating system.");
        }

        try ( final InputStream is = MiniMap2Index.class.getResourceAsStream(libName) ) {
            if ( is == null ) {
                throw new MiniMap2Exception("Can't find resource " + libName);
            }
            final File tmpFile = File.createTempFile("libbwa.",".jnilib");
            tmpFile.deleteOnExit();
            Files.copy(is, tmpFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
            System.load(tmpFile.getPath());
        } catch ( final IOException ioe ) {
            throw new MiniMap2Exception("Misconfiguration: Unable to load minimap2 native library " + libName, ioe);
        }
        nativeLibLoaded = true;
    }

    /* returns a ByteBuffer containing this structure:
        typedef struct {
        int64_t flag;    // see MM_F_* macros
        int seed;
        int sdust_thres; // score threshold for SDUST; 0 to disable

        int max_qlen;    // max query length

        int bw;          // bandwidth
        int max_gap, max_gap_ref; // break a chain if there are no minimizers in a max_gap window
        int max_frag_len;
        int max_chain_skip, max_chain_iter;
        int min_cnt;         // min number of minimizers on each chain
        int min_chain_score; // min chaining score

        float mask_level;
        float pri_ratio;
        int best_n;      // top best_n chains are subjected to DP alignment

        int max_join_long, max_join_short;
        int min_join_flank_sc;
        float min_join_flank_ratio;

        int a, b, q, e, q2, e2; // matching score, mismatch, gap-open and gap-ext penalties
        int sc_ambi; // score when one or both bases are "N"
        int noncan;      // cost of non-canonical splicing sites
        int junc_bonus;
        int zdrop, zdrop_inv;   // break alignment if alignment score drops too fast along the diagonal
        int end_bonus;
        int min_dp_max;  // drop an alignment if the score of the max scoring segment is below this threshold
        int min_ksw_len;
        int anchor_ext_len, anchor_ext_shift;
        float max_clip_ratio; // drop an alignment if BOTH ends are clipped above this ratio

        int pe_ori, pe_bonus;

        float mid_occ_frac;  // only used by mm_mapopt_update(); see below
        int32_t min_mid_occ;
        int32_t mid_occ;     // ignore seeds with occurrences above this threshold
        int32_t max_occ;
        int mini_batch_size; // size of a batch of query bases to process in parallel
        int64_t max_sw_mat;

        const char *split_prefix;
      } mm_mapopt_t;
    */
    static native ByteBuffer createOptions( long addr, String preset );

    // allocates native memory and returns it in a ByteBuffer
    static native ByteBuffer createSeqBuffer( int length );

    // frees native memory
    static native void destroyByteBuffer( ByteBuffer buffer );

    public static native String getVersion();

    private static native long openIndex( String mmiFile );

    private static native void destroyIndex( long addr );

    /*
     addr is the native memory address of the index
     opts is a ByteBuffer than contains a mm_mapopt_t structure
     seqs is a ByteBuffer that contains:
       a 32-bit integer count of the number of sequences to follow
       a 32-bit integer for each sequence, giving its length
       the base calls for each sequence, all run together in a single pool of bytes

     we return a ByteBuffer that contains:
     for each sequence,
       a 32-bit integer count of the number of alignments that follow
       for each alignment, a pseudo-structure like this:
      typedef struct {
          int32_t refID; // reference id
          int32_t pos;   // reference starting position (0-based) -- if negative, it means reverse strand starting at ~pos
          int32_t mapQ;  // the map quality
          int32_t nCigar; // nCigarOps
          int32_t cigarOp[nCigarOps]; // len<<4 | op (i.e., the usual BAM encoding)
      } Alignment;
    */
    private static native ByteBuffer createAlignments( long addr, ByteBuffer opts, ByteBuffer seqs );

    // returns a ByteBuffer with all the reference contig names concatenated (null byte delimited)
    private static native ByteBuffer getRefNames( long addr );
}
