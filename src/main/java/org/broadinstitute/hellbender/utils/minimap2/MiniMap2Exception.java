package org.broadinstitute.hellbender.utils.minimap2;

public class MiniMap2Exception extends RuntimeException {
    public MiniMap2Exception( final String message ) { super(message); }
    public MiniMap2Exception( final String message, final Exception cause ) { super(message, cause); }
}
