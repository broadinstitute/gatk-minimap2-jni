package org.broadinstitute.hellbender.utils.minimap2;

import java.nio.ByteBuffer;
import java.util.Objects;

public class MiniMap2Alignment {
    private final int refId;
    private final int refStart;
    private final boolean revStrand;
    private final int mapQ;
    private final String cigar;

    public MiniMap2Alignment( final ByteBuffer alignBuffer ) {
        refId = alignBuffer.getInt();
        refStart = alignBuffer.getInt();
        final int mapQAndStrand = alignBuffer.getInt();
        if ( mapQAndStrand < 0 ) {
            revStrand = true;
            mapQ = ~mapQAndStrand;
        } else {
            revStrand = false;
            mapQ = mapQAndStrand;
        }

        final StringBuilder sb = new StringBuilder();
        int nCigarWords = alignBuffer.getInt();
        while ( nCigarWords-- > 0 ) {
            final int cigarWord = alignBuffer.getInt();
            final char operator = "MIDNSHP=X".charAt(cigarWord & 0x0f);
            sb.append(cigarWord >>> 4).append(operator);
        }
        cigar = sb.toString();
    }

    public MiniMap2Alignment( final int refId, final int refStart, final boolean revStrand,
                              final int mapQ, final String cigar ) {
        this.refId = refId;
        this.refStart = refStart;
        this.revStrand = revStrand;
        this.mapQ = mapQ;
        this.cigar = cigar;
    }

    public int getRefId() { return refId; }
    public int getRefStart() { return refStart; }
    public boolean isRevStrand() { return revStrand; }
    public int getMapQ() { return mapQ; }
    public String getCigar() { return cigar; }

    @Override public boolean equals( final Object obj ) {
        if ( obj == this ) return true;
        if ( !(obj instanceof MiniMap2Alignment) ) return false;
        final MiniMap2Alignment that = (MiniMap2Alignment)obj;
        return this.refId == that.refId && this.refStart == that.refStart && this.revStrand == that.revStrand &&
                this.mapQ == that.mapQ && Objects.equals(this.cigar, that.cigar);
    }
}
