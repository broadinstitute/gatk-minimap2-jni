#include <jni.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include "version.h"
#include "minimap2/minimap.h"

/*
 * Implementation of native routines declared in MiniMap2Index.java.
    static native ByteBuffer createOptions( String preset );
    static native ByteBuffer createSeqBuffer( int length );
    static native void destroyByteBuffer( ByteBuffer buffer );

    private static native long openIndex( String mmiFile );
    private static native void destroyIndex( long indexAddress );
    private static native ByteBuffer align( long addr, ByteBuffer seqs, ByteBuffer opts );
    private static native ByteBuffer getRefNames( long addr );

    private static native String getVersion();
 */

static char* jstringToChars( JNIEnv* env, jstring in ) {
    const char* tmp = (*env)->GetStringUTFChars(env, in, 0);
    char* res = strdup(tmp);
    (*env)->ReleaseStringUTFChars(env, in, tmp);
    return res;
}

static void throwErrorMessage( JNIEnv* env, char* message ) {
    jclass exceptionClass = (*env)->FindClass(env, "org/broadinstitute/hellbender/utils/minimap2/MiniMap2Exception");
    if ( exceptionClass ) (*env)->ThrowNew(env, exceptionClass, message);
}

JNIEXPORT jobject JNICALL
Java_org_broadinstitute_hellbender_utils_minimap2_MiniMap2Index_createOptions( JNIEnv* env, jclass cls, jlong idxAddr, jstring preset ) {
    mm_idxopt_t idxOpts;
    mm_mapopt_t* pMapOpts = malloc(sizeof(mm_mapopt_t));
    if ( !pMapOpts ) {
        throwErrorMessage(env, "C code can't allocate memory for options buffer");
        return 0;
    }
    mm_set_opt(0, &idxOpts, pMapOpts);
    if ( preset ) {
        char* presetChars = jstringToChars(env, preset);
        mm_set_opt(presetChars, &idxOpts, pMapOpts);
        free(presetChars);
    }
    pMapOpts->flag |= MM_F_CIGAR;
    mm_mapopt_update(pMapOpts, (mm_idx_t*)idxAddr);
    return (*env)->NewDirectByteBuffer(env, pMapOpts, sizeof(mm_mapopt_t));
}

JNIEXPORT jobject JNICALL
Java_org_broadinstitute_hellbender_utils_minimap2_MiniMap2Index_createSeqBuffer( JNIEnv* env, jclass cls, jint length ) {
    void* buf = malloc(length);
    if ( !buf ) {
        throwErrorMessage(env, "C code can't allocate memory for sequence buffer");
        return 0;
    }
    return (*env)->NewDirectByteBuffer(env, buf, length);
}

JNIEXPORT void JNICALL
Java_org_broadinstitute_hellbender_utils_minimap2_MiniMap2Index_destroyByteBuffer( JNIEnv* env, jclass cls, jobject byteBuf ) {
    void* buf = (*env)->GetDirectBufferAddress(env, byteBuf);
    if ( !buf ) throwErrorMessage(env, "C code can't get ByteBuffer address");
    free(buf);
}

JNIEXPORT jlong JNICALL
Java_org_broadinstitute_hellbender_utils_minimap2_MiniMap2Index_openIndex( JNIEnv* env, jclass cls, jstring idxFilename ) {
    mm_verbose = 0;
    char* imageName = jstringToChars(env, idxFilename);
    mm_idx_reader_t* pIdxReader = mm_idx_reader_open(imageName, 0, 0);
    free(imageName);
    if ( !pIdxReader ) {
        throwErrorMessage(env, "C code can't open index file");
        return 0;
    }
    mm_idx_t* pIdx = mm_idx_reader_read(pIdxReader, 0);
    if ( !pIdx ) throwErrorMessage(env, "C code can't read index file");
    mm_idx_reader_close(pIdxReader);
    return (jlong)pIdx;
}

JNIEXPORT void JNICALL
Java_org_broadinstitute_hellbender_utils_minimap2_MiniMap2Index_destroyIndex( JNIEnv* env, jclass cls, jlong idxAddr ) {
    if ( !idxAddr ) {
        throwErrorMessage(env, "C code can't close null index address");
        return;
    }
    mm_idx_destroy((mm_idx_t*)idxAddr);
}

// we accept a ByteBuffer that contains:
//   a 32-bit integer count of the number of sequences to follow
//   a 32-bit integer for each sequence, giving its length
//   the base calls for each sequence, all run together in a single pool of bytes
// we return a ByteBuffer that contains:
// for each sequence,
//   a 32-bit integer count of the number of alignments that follow
//   for each alignment, a pseudo-structure like this:
/*
typedef struct {
    int32_t refID; // reference id
    int32_t pos;   // reference starting position (0-based) -- if negative, it means reverse strand starting at ~pos
    int32_t mapQ;  // the map quality
    int32_t nCigar; // nCigarOps
    int32_t cigarOp[nCigarOps]; // len<<4 | op (i.e., the usual BAM encoding)
} Alignment;
*/
JNIEXPORT jobject JNICALL
Java_org_broadinstitute_hellbender_utils_minimap2_MiniMap2Index_createAlignments(
                JNIEnv* env, jclass cls, jlong idxAddr, jobject optsBuf, jobject seqsBuf ) {
    if ( !idxAddr ) {
        throwErrorMessage(env, "C code can't align with a null index address");
        return 0;
    }
    mm_idx_t* pIdx = (mm_idx_t*)idxAddr;
    mm_mapopt_t* pOpts = (*env)->GetDirectBufferAddress(env, optsBuf);
    if ( !pOpts ) {
        throwErrorMessage(env, "C code can't get address for opts ByteBuffer");
        return 0;
    }
    uint32_t* pLengths = (*env)->GetDirectBufferAddress(env, seqsBuf);
    if ( !pLengths ) {
        throwErrorMessage(env, "C code can't get address for seqs ByteBuffer");
        return 0;
    }

    uint32_t nSeqs = *pLengths++;
    char* pSeqs = (char*)(pLengths + nSeqs);
    mm_tbuf_t* pTBuf = mm_tbuf_init();
    uint32_t** allAlignsBase = calloc(nSeqs, sizeof(uint32_t*));
    if ( !allAlignsBase ) {
        throwErrorMessage(env, "C code can't allocate temporary memory for alignments");
        return 0;
    }
    uint32_t** ppAligns = allAlignsBase;
    uint32_t** ppAlignsEnd = ppAligns + nSeqs;
    int seqId;
    int nAligns;
    for ( seqId = 0; seqId != nSeqs; ++seqId ) {
        uint32_t seqLen = *pLengths++;
        nAligns = 0;
        mm_reg1_t* pAlignsBase = mm_map(pIdx, seqLen, pSeqs, &nAligns, pTBuf, pOpts, 0);
        pSeqs += seqLen;
        size_t len = sizeof(uint32_t); // space for nAligns
        mm_reg1_t* pAlign = pAlignsBase;
        mm_reg1_t* pEnd = pAlign + nAligns;
        while ( pAlign != pEnd ) {
            uint32_t nCigar = pAlign->p ? (pAlign->p->n_cigar + (pAlign->qs > 0) + (pAlign->qe < seqLen)) : 0;
            len += (4 + nCigar) * sizeof(uint32_t);
            pAlign += 1;
        }
        uint32_t* buf = malloc(len + sizeof(uint32_t)); // extra space for buffer len
        *ppAligns++ = buf;
        if ( buf ) {
            *buf++ = len;
            *buf++ = nAligns;
        }
        pAlign = pAlignsBase;
        while ( pAlign != pEnd ) {
            if ( buf ) {
                *buf++ = pAlign->rid;
                *buf++ = pAlign->rs;
                *buf++ = pAlign->rev ? ~pAlign->mapq : pAlign->mapq;
                *buf++ = pAlign->p ? (pAlign->p->n_cigar + (pAlign->qs > 0) + (pAlign->qe < seqLen)) : 0;
                if ( pAlign->qs > 0 ) {
                    *buf++ = (pAlign->qs << 4) | 4; // soft-clip at the beginning
                }
                if ( pAlign->p ) {
                    uint32_t nCigar = pAlign->p->n_cigar;
                    memcpy(buf, pAlign->p->cigar, nCigar * sizeof(uint32_t));
                    buf += nCigar;
                }
                if ( pAlign->qe < seqLen ) {
                    *buf++ = ((seqLen - pAlign->qe) << 4) | 4; // final soft-clip
                }
            }
            free(pAlign->p);
            pAlign += 1;
        }
        free(pAlignsBase);
        if ( !buf ) break;
    }
    mm_tbuf_destroy(pTBuf);
    ppAligns = allAlignsBase;
    uint32_t totLen = 0;
    while ( ppAligns != ppAlignsEnd ) {
        uint32_t* pBuf = *ppAligns++;
        if ( !pBuf ) {
            totLen = 0;
            break;
        }
        totLen += *pBuf;
    }
    char* allBufBase = totLen ? malloc(totLen) : 0;
    char* allBuf = allBufBase;
    ppAligns = allAlignsBase;
    while ( ppAligns != ppAlignsEnd ) {
        uint32_t* pBuf = *ppAligns++;
        int len = *pBuf;
        if ( allBuf ) {
            memcpy(allBuf, pBuf + 1, len);
            allBuf += len;
        }
        free(pBuf);
    }
    ppAligns = ppAlignsEnd;
    while ( ppAligns != ppAlignsEnd ) {
        free(*ppAligns++);
    }
    free(allAlignsBase);
    if ( !allBufBase ) {
        throwErrorMessage(env, "C code can't create memory for alignment buffer");
        return 0;
    }
    jobject alnBuf = (*env)->NewDirectByteBuffer(env, allBufBase, totLen);
    if ( !alnBuf ) {
        free(allBuf);
        throwErrorMessage(env, "C code can't create ByteBuffer for alignments");
    }
    return alnBuf;
}

// returns a ByteBuffer with the reference contig names concatenated
JNIEXPORT jobject JNICALL
Java_org_broadinstitute_hellbender_utils_minimap2_MiniMap2Index_getRefNames( JNIEnv* env, jclass cls, jlong idxAddr ) {
    if ( !idxAddr ) {
        throwErrorMessage(env, "C code can't get ref names from null index address");
        return 0;
    }
    mm_idx_t* pIdx = (mm_idx_t*)idxAddr;
    uint32_t nSeqs = pIdx->n_seq;
    mm_idx_seq_t* pSeqs = pIdx->seq;
    mm_idx_seq_t* pEnd = pSeqs + nSeqs;
    uint32_t bufLen = nSeqs; // a null byte for each name
    while ( pSeqs != pEnd ) bufLen += strlen(pSeqs++->name);
    char* bufBase = malloc(bufLen);
    if ( !bufBase ) {
        throwErrorMessage(env, "C code can't allocate memory for ref names");
        return 0;
    }
    char* buf = bufBase;
    pSeqs = pIdx->seq;
    while ( pSeqs != pEnd ) {
        char* name = pSeqs++->name;
        uint32_t len = strlen(name) + 1;
        memcpy(buf, name, len);
        buf += len;
    }
    jobject namesBuf = (*env)->NewDirectByteBuffer(env, bufBase, bufLen);
    if ( !namesBuf ) {
        free(bufBase);
        throwErrorMessage(env, "C code unable to create ByteBuffer for ref names");
    }
    return namesBuf;
}

JNIEXPORT jstring JNICALL
Java_org_broadinstitute_hellbender_utils_minimap2_MiniMap2Index_getVersion( JNIEnv* env, jclass cls ) {
    return (*env)->NewStringUTF(env, MINIMAP2_VERSION);
}
