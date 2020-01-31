package org.broadinstitute.hellbender.utils.minimap2;

import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.nio.ByteBuffer;
import java.util.Arrays;

@Test
public class MiniMap2IndexUnitTest {
    private static MiniMap2Index index;

    @BeforeClass
    void openIndex() {
        index = new MiniMap2Index("src/test/resources/org/broadinstitute/hellbender/utils/minimap2/test.mmi");
    }

    @AfterClass
    void closeIndex() {
        index.close(); index = null;
    }

    @Test
    void testOpen() {
        Assert.assertTrue(index.isOpen());
    }

    @Test
    void testOptsSize() {
        final ByteBuffer optsBuf = index.createOptions(MiniMap2Aligner.Preset.ASM5.getName());
        Assert.assertEquals(optsBuf.capacity(), MiniMap2Aligner.EXPECTED_OPTS_SIZE);
        MiniMap2Index.destroyByteBuffer(optsBuf);
    }

    @Test
    void testGetRefNames() {
        Assert.assertEquals(index.getRefNames(), Arrays.asList("ref1", "ref2"));
    }

    @Test
    void testGetVersion() {
        final String version = MiniMap2Index.getVersion();
        Assert.assertTrue( version != null && !version.isEmpty());
    }
}
