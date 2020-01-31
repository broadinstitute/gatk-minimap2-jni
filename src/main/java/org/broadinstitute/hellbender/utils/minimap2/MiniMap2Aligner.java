package org.broadinstitute.hellbender.utils.minimap2;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.function.Function;

/**
 * Given an open index, this class lets you do alignment of sequences.
 * Don't forget to close it, or you'll leak a little memory.
 * Usage pattern:
 *   Create a MiniMap2Aligner on some MiniMap2Index, specifying a preset if desired
 *   Tweak options as necessary using get/set methods.
 *   Align 1 or more chunks of sequences with alignSeqs
 *   Close the MiniMap2Aligner
 * This class is not thread-safe, but it's very light-weight:  just use a separate instance in each thread.
 */
public class MiniMap2Aligner implements AutoCloseable {
	private final MiniMap2Index index;
	private ByteBuffer opts;
	static final int EXPECTED_OPTS_SIZE = 192;

	public MiniMap2Aligner( final MiniMap2Index index ) { this(index, null); }

	public MiniMap2Aligner( final MiniMap2Index index, final Preset preset ) {
		if ( !index.isOpen() ) {
			throw new MiniMap2Exception("Can't create MiniMap2Aligner: The index has been closed.");
		}
		this.index = index;
		opts = index.createOptions(preset == null ? null : preset.getName());
		if ( opts == null ) {
			throw new MiniMap2Exception("Can't create MiniMap2Aligner: Unable to retrieve options object.");
		}
		opts.order(ByteOrder.nativeOrder()).clear();
		if ( opts.capacity() != EXPECTED_OPTS_SIZE ) {
			close();
			throw new MiniMap2Exception(
					"Can't create MiniMap2Aligner: Unexpected options object size suggests wrong version of MiniMap2.");
		}
	}

	/**
	 * Align some sequences.
	 * @param sequences A list of byte[]'s that contain base calls (ASCII 'A', 'C', 'G', or 'T').
	 * @return A list of the same length as the input list.  Each element is a list of alignments for the corresponding sequence.
	 */
	public List<List<MiniMap2Alignment>> alignSeqs( final List<byte[]> sequences ) {
		return alignSeqs(sequences, seq -> seq);
	}

	/**
	 * A more abstract version that takes an iterable of things that can be turned into a byte[] of base calls.
	 * @param iterable An iterable over something like a read, that contains a sequence.
	 * @param func A lambda that picks the sequence out of your read-like thing.
	 * @param <T> The read-like thing.
	 * @return A list of (possibly multiple) alignments for each input sequence.
	 */
	public <T> List<List<MiniMap2Alignment>> alignSeqs( final Iterable<T> iterable, final Function<T, byte[]> func ) {
		int nSeqs = 0;
		int nBases = 0;
		for ( final T item : iterable ) {
			nSeqs += 1;
			nBases += func.apply(item).length;
		}

		final int len = 4 + nSeqs * 4 + nBases;
		final ByteBuffer seqBuffer = MiniMap2Index.createSeqBuffer(len);
		if ( seqBuffer == null ) {
			throw new MiniMap2Exception("Can't create buffer for passing sequences to MiniMap2.");
		}
		seqBuffer.order(ByteOrder.nativeOrder()).clear();

		seqBuffer.putInt(nSeqs);
		for ( final T item : iterable ) {
			seqBuffer.putInt(func.apply(item).length);
		}
		for ( final T item : iterable ) {
			seqBuffer.put(func.apply(item));
		}

		try {
			final ByteBuffer alignBuffer = index.align(opts, seqBuffer);
			if ( alignBuffer == null ) {
				throw new MiniMap2Exception("Couldn't create alignments.  Reasons are obscure.");
			}
			try {
				alignBuffer.order(ByteOrder.nativeOrder()).position(0).limit(alignBuffer.capacity());
				final List<List<MiniMap2Alignment>> result = new ArrayList<>(nSeqs);
				while ( nSeqs-- > 0 ) {
					int nAligns = alignBuffer.getInt();
					final List<MiniMap2Alignment> aligns = new ArrayList<>(nAligns);
					while ( nAligns-- > 0 ) {
						aligns.add(new MiniMap2Alignment(alignBuffer));
					}
					result.add(aligns);
				}
				return result;
			} finally {
				MiniMap2Index.destroyByteBuffer(alignBuffer);
			}
		} finally {
			MiniMap2Index.destroyByteBuffer(seqBuffer);
		}
	}

	public boolean isOpen() { return opts != null; }

	public MiniMap2Index getIndex() { return index; }

	@Override
	public void close() {
		synchronized (this) {
			if ( opts != null ) {
				MiniMap2Index.destroyByteBuffer(opts);
				opts = null;
			}
		}
	}

	/** sets of options for various types of data */
	public enum Preset {
		AVA_ONT("ava-ont"),
		AVA_PB("ava-pb"),
		MAP_10K("map10k"),
		MAP_PB("map-pb"),
		MAP_ONT("map-ont"),
		ASM5("asm5"),
		ASM10("asm10"),
		ASM20("asm20"),
		SHORT("short"),
		SR("sr"),
		SPLICE("splice"),
		CDNA("cdna");

		private final String name;
		Preset( final String name ) { this.name = name; }
		public String getName() { return name; }
	}

	// bit-wise constants for get/setFlag methods
	public static final long MM_F_NO_DIAG       = 0x001; // no exact diagonal hit
	public static final long MM_F_NO_DUAL       = 0x002; // skip pairs where query name is lexicographically larger than target name
	public static final long MM_F_CIGAR         = 0x004;
	public static final long MM_F_OUT_SAM       = 0x008;
	public static final long MM_F_NO_QUAL       = 0x010;
	public static final long MM_F_OUT_CG        = 0x020;
	public static final long MM_F_OUT_CS        = 0x040;
	public static final long MM_F_SPLICE        = 0x080; // splice mode
	public static final long MM_F_SPLICE_FOR    = 0x100; // match GT-AG
	public static final long MM_F_SPLICE_REV    = 0x200; // match CT-AC, the reverse complement of GT-AG
	public static final long MM_F_NO_LJOIN      = 0x400;
	public static final long MM_F_OUT_CS_LONG   = 0x800;
	public static final long MM_F_SR            = 0x1000;
	public static final long MM_F_FRAG_MODE     = 0x2000;
	public static final long MM_F_NO_PRINT_2ND  = 0x4000;
	public static final long MM_F_2_IO_THREADS  = 0x8000;
	public static final long MM_F_LONG_CIGAR    = 0x10000;
	public static final long MM_F_INDEPEND_SEG  = 0x20000;
	public static final long MM_F_SPLICE_FLANK  = 0x40000;
	public static final long MM_F_SOFTCLIP      = 0x80000;
	public static final long MM_F_FOR_ONLY      = 0x100000;
	public static final long MM_F_REV_ONLY      = 0x200000;
	public static final long MM_F_HEAP_SORT     = 0x400000;
	public static final long MM_F_ALL_CHAINS    = 0x800000;
	public static final long MM_F_OUT_MD        = 0x1000000;
	public static final long MM_F_COPY_COMMENT  = 0x2000000;
	public static final long MM_F_EQX           = 0x4000000; // use =/X instead of M
	public static final long MM_F_PAF_NO_HIT    = 0x8000000; // output unmapped reads to PAF
	public static final long MM_F_NO_END_FLT    = 0x10000000;
	public static final long MM_F_HARD_MLEVEL   = 0x20000000;
	public static final long MM_F_SAM_HIT_ONLY  = 0x40000000;

	public long getFlags() { return getOpts().getLong(0); }
	public void setFlags( final long arg ) { getOpts().putLong(0, arg); }

	public int getSeed() { return getOpts().getInt(8); }
	public void setSeed( final int arg ) { getOpts().putInt(8, arg); }
	public int getSDustThreshold() { return getOpts().getInt(12); }
	public void setSDustThreshold( final int arg ) { getOpts().putInt(12, arg); }
	public int getMaxQueryLen() { return getOpts().getInt(16); }
	public void setMaxQueryLen( final int arg ) { getOpts().putInt(16, arg); }
	public int getBandwidth() { return getOpts().getInt(20); }
	public void setBandwidth( final int arg ) { getOpts().putInt(20, arg); }
	// break a chain if there are no minimizers in a max_gap window
	public int getMaxGap() { return getOpts().getInt(24); }
	public void setMaxGap( final int arg ) { getOpts().putInt(24, arg); }
	public int getMaxGapRef() { return getOpts().getInt(28); }
	public void setMaxGapRef( final int arg ) { getOpts().putInt(28, arg); }
	public int getMaxFragLen() { return getOpts().getInt(32); }
	public void setMaxFragLen( final int arg ) { getOpts().putInt(32, arg); }
	public int getMaxChainSkip() { return getOpts().getInt(36); }
	public void setMaxChainSkip( final int arg ) { getOpts().putInt(36, arg); }
	public int getMaxChainIter() { return getOpts().getInt(40); }
	public void setMaxChainIter( final int arg ) { getOpts().putInt(40, arg); }
	// min number of minimizers on each chain
	public int getMinCnt() { return getOpts().getInt(44); }
	public void setMinCnt( final int arg ) { getOpts().putInt(44, arg); }
	public int getMinChainScore() { return getOpts().getInt(48); }
	public void setMinChainScore( final int arg ) { getOpts().putInt(48, arg); }
	public float getMaskLevel() { return getOpts().getFloat(52); }
	public void setMaskLevel( final float arg ) { getOpts().putFloat(52, arg); }
	public float getPriRatio() { return getOpts().getFloat(56); }
	public void setPriRatio( final float arg ) { getOpts().putFloat(56, arg); }
	// top best_n chains are subjected to DP alignment
	public int getBestN() { return getOpts().getInt(60); }
	public void setBestN( final int arg ) { getOpts().putInt(60, arg); }
	public int getMaxJoinLong() { return getOpts().getInt(64); }
	public void setMaxJoinLong( final int arg ) { getOpts().putInt(64, arg); }
	public int getMaxJoinShort() { return getOpts().getInt(68); }
	public void setMaxJoinShort( final int arg ) { getOpts().putInt(68, arg); }
	public int getMinJoinFlankScore() { return getOpts().getInt(72); }
	public void setMinJoinFlankScore( final int arg ) { getOpts().putInt(72, arg); }
	public float getMinJoinFlankRatio() { return getOpts().getFloat(76); }
	public void setMinJoinFlankRatio( final float arg ) { getOpts().putFloat(76, arg); }
	// matching score
	public int getA() { return getOpts().getInt(80); }
	public void setA( final int arg ) { getOpts().putInt(80, arg); }
	// mismatch score
	public int getB() { return getOpts().getInt(84); }
	public void setB( final int arg ) { getOpts().putInt(84, arg); }
	// gap open penalty
	public int getQ() { return getOpts().getInt(88); }
	public void setQ( final int arg ) { getOpts().putInt(88, arg); }
	// gap extend penalty
	public int getE() { return getOpts().getInt(92); }
	public void setE( final int arg ) { getOpts().putInt(92, arg); }
	public int getQ2() { return getOpts().getInt(96); }
	public void setQ2( final int arg ) { getOpts().putInt(96, arg); }
	public int getE2() { return getOpts().getInt(100); }
	public void setE2( final int arg ) { getOpts().putInt(100, arg); }
	// score when one or both bases are "N"
	public int getScoreAmbi() { return getOpts().getInt(104); }
	public void setScoreAmbi( final int arg ) { getOpts().putInt(104, arg); }
	// cost of non-canonical splicing sites
	public int getNonCan() { return getOpts().getInt(108); }
	public void setNonCan( final int arg ) { getOpts().putInt(108, arg); }
	public int getJuncBonus() { return getOpts().getInt(112); }
	public void setJuncBonus( final int arg ) { getOpts().putInt(112, arg); }
	// break alignment if alignment score drops too fast along the diagonal
	public int getZDrop() { return getOpts().getInt(116); }
	public void setZDrop( final int arg ) { getOpts().putInt(116, arg); }
	public int getZDropInv() { return getOpts().getInt(120); }
	public void setZDropInv( final int arg ) { getOpts().putInt(120, arg); }
	public int getEndBonus() { return getOpts().getInt(124); }
	public void setEndBonus( final int arg ) { getOpts().putInt(124, arg); }
	// drop an alignment if the score of the max scoring segment is below this threshold
	public int getMinDPMax() { return getOpts().getInt(128); }
	public void setMinDPMax( final int arg ) { getOpts().putInt(128, arg); }
	public int getMinKSWLen() { return getOpts().getInt(132); }
	public void setMinKSWLen( final int arg ) { getOpts().putInt(132, arg); }
	public int getAnchorExtLen() { return getOpts().getInt(136); }
	public void setAnchorExtLen( final int arg ) { getOpts().putInt(136, arg); }
	public int getAnchorExtShift() { return getOpts().getInt(140); }
	public void setAnchorExtShift( final int arg ) { getOpts().putInt(140, arg); }
	// drop an alignment if BOTH ends are clipped above this ratio
	public float getMaxClipRatio() { return getOpts().getFloat(144); }
	public void setMaxClipRatio( final float arg ) { getOpts().putFloat(144, arg); }
	public int getPEOri() { return getOpts().getInt(148); }
	public void setPEOri( final int arg ) { getOpts().putInt(148, arg); }
	public int getPEBonus() { return getOpts().getInt(152); }
	public void setPEBonus( final int arg ) { getOpts().putInt(152, arg); }
	public float getMidOccFrac() { return getOpts().getFloat(156); }
	public void setMidOccFrac( final float arg ) { getOpts().putFloat(156, arg); }
	public int getMinMidOcc() { return getOpts().getInt(160); }
	public void setMinMidOcc( final int arg ) { getOpts().putInt(160, arg); }
	// ignore seeds with occurrences above this threshold
	public int getMidOcc() { return getOpts().getInt(164); }
	public void setMidOcc( final int arg ) { getOpts().putInt(164, arg); }
	public int getMaxOcc() { return getOpts().getInt(168); }
	public void setMaxOcc( final int arg ) { getOpts().putInt(168, arg); }
	// size of a batch of query bases to process in parallel
	public int getMiniBatchSize() { return getOpts().getInt(172); }
	public void setMiniBatchSize( final int arg ) { getOpts().putInt(172, arg); }
	public long getMaxSWMat() { return getOpts().getLong(176); }
	public void setMaxSWMat( final long arg ) { getOpts().putLong(176, arg); }
	// split_prefix is ignored

	public String toString() {
		return "Flags: " + getOpts().getLong(0) +
			"\nSeed: " + getOpts().getInt(8) +
			"\nSDustThreshold: " + getOpts().getInt(12) +
		    "\nMaxQueryLen: " + getOpts().getInt(16) +
		    "\nBandwidth: " + getOpts().getInt(20) +
		    "\nMaxGap: " + getOpts().getInt(24) +
		    "\nMaxGapRef: " + getOpts().getInt(28) +
		    "\nMaxFragLen: " + getOpts().getInt(32) +
		    "\nMaxChainSkip: " + getOpts().getInt(36) +
		    "\nMaxChainIter: " + getOpts().getInt(40) +
		    "\nMinCnt: " + getOpts().getInt(44) +
		    "\nMinChainScore: " + getOpts().getInt(48) +
		    "\nMaskLevel: " + getOpts().getFloat(52) +
		    "\nPriRatio: " + getOpts().getFloat(56) +
		    "\nBestN: " + getOpts().getInt(60) +
		    "\nMaxJoinLong: " + getOpts().getInt(64) +
		    "\nMaxJoinShort: " + getOpts().getInt(68) +
		    "\nMinJoinFlankScore: " + getOpts().getInt(72) +
		    "\nMinJoinFlankRatio: " + getOpts().getFloat(76) +
		    "\nA: " + getOpts().getInt(80) +
		    "\nB: " + getOpts().getInt(84) +
		    "\nQ: " + getOpts().getInt(88) +
		    "\nE: " + getOpts().getInt(92) +
		    "\nQ2: " + getOpts().getInt(96) +
		    "\nE2: " + getOpts().getInt(100) +
		    "\nScoreAmbi: " + getOpts().getInt(104) +
		    "\nNonCan: " + getOpts().getInt(108) +
		    "\nJuncBonus: " + getOpts().getInt(112) +
		    "\nZDrop: " + getOpts().getInt(116) +
		    "\nZDropInv: " + getOpts().getInt(120) +
		    "\nEndBonus: " + getOpts().getInt(124) +
		    "\nMinDPMax: " + getOpts().getInt(128) +
		    "\nMinKSWLen: " + getOpts().getInt(132) +
		    "\nAnchorExtLen: " + getOpts().getInt(136) +
		    "\nAnchorExtShift: " + getOpts().getInt(140) +
		    "\nMaxClipRatio: " + getOpts().getFloat(144) +
		    "\nPEOri: " + getOpts().getInt(148) +
		    "\nPEBonus: " + getOpts().getInt(152) +
		    "\nMidOccFrac: " + getOpts().getFloat(156) +
		    "\nMinMidOcc: " + getOpts().getInt(160) +
		    "\nMidOcc: " + getOpts().getInt(164) +
		    "\nMaxOcc: " + getOpts().getInt(168) +
		    "\nMiniBatchSize: " + getOpts().getInt(172) +
		    "\nMaxSWMat: " + getOpts().getLong(176);
	}
	private ByteBuffer getOpts() {
		if ( opts == null ) {
			throw new IllegalStateException("The aligner has been closed.");
		}
		return opts;
	}
}
