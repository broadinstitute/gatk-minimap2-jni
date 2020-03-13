package org.broadinstitute.hellbender.utils.minimap2;

import java.nio.ByteBuffer;
import java.util.Objects;

public class MiniMap2Alignment {
    private final int samFlag;
    private final int refId;
    private final int refStart;
    private final int mapQ;
    private final int nm;
    private final String cigar;

    public MiniMap2Alignment( final ByteBuffer alignBuffer ) {
        samFlag = alignBuffer.getInt();
        refId = alignBuffer.getInt();
        refStart = alignBuffer.getInt();
        mapQ = alignBuffer.getInt();

        final StringBuilder sb = new StringBuilder();
        int nCigarWords = alignBuffer.getInt();
        nm = nCigarWords > 0 ? alignBuffer.getInt() : 0;
        while ( nCigarWords-- > 0 ) {
            final int cigarWord = alignBuffer.getInt();
            final char operator = "MIDNSHP=X".charAt(cigarWord & 0x0f);
            sb.append(cigarWord >>> 4).append(operator);
        }
        cigar = sb.toString();
    }

    public MiniMap2Alignment( final int samFlag, final int refId, final int refStart,
                              final int mapQ, final int nm, final String cigar ) {
        this.samFlag = samFlag;
        this.refId = refId;
        this.refStart = refStart;
        this.mapQ = mapQ;
        this.nm = nm;
        this.cigar = cigar;
    }

    public int getSAMFlag() { return samFlag; }
    public int getRefId() { return refId; }
    public int getRefStart() { return refStart; }
    public int getMapQ() { return mapQ; }
    public int getNM() { return nm; }
    public String getCigar() { return cigar; }

    @Override public boolean equals( final Object obj ) {
        if ( obj == this ) return true;
        if ( !(obj instanceof MiniMap2Alignment) ) return false;
        final MiniMap2Alignment that = (MiniMap2Alignment)obj;
        return this.samFlag == that.samFlag && this.refId == that.refId && this.refStart == that.refStart &&
                this.mapQ == that.mapQ && this.nm == that.nm && Objects.equals(this.cigar, that.cigar);
    }
}
