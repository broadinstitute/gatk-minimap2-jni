package org.broadinstitute.hellbender.utils.minimap2;

import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

@Test
public class MiniMap2AlignerUnitTest {
    private static MiniMap2Index index;

    @BeforeClass
    void openIndex() {
        index = new MiniMap2Index("src/test/resources/org/broadinstitute/hellbender/utils/minimap2/test.mmi");
    }

    @AfterClass
    void closeIndex() { index.close(); index = null; }

    @Test
    void testAlignment() {
        try ( final MiniMap2Aligner aligner = new MiniMap2Aligner(index) ) {
            final List<List<MiniMap2Alignment>> alignments =
                    aligner.alignSeqs(Arrays.asList(
                    "AGAACTCCACACATGGGATAATGTTTTGGCTAGATGGCTCCCCTACTTAGAAACATACAATTGCTAGTCATATTTCTAATTTTAGGATTTCGAGATACTGGTGATGAAGATCACATGTCC".getBytes(),
                    "TTTTTTTTACACATGGGATAATGTTTTGGCTAGATGGCTCCCCTACTTAGAAACATACAATTGCTAGTCATATTTCTAATTTTAGGATTTCGAGATACTGGTGATGAAGATCACATGTCC".getBytes(),
                    "AGAACTCCACACATGGGATAATGTTTTGGCTAGATGGCTCCCCTACTTAGAAACATACAATTGCTAGTCATATTTCTAATTTTAGGATTTCGAGATACTGGTGATGAAGATCGGGGGGGG".getBytes(),
                    "TTTTTTTTACACATGGGATAATGTTTTGGCTAGATGGCTCCCCTACTTAGAAACATACAATTGCTAGTCATATTTCTAATTTTAGGATTTCGAGATACTGGTGATGAAGATCGGGGGGGG".getBytes(),
                    "ATCCAAAGAAAAGCAGAGAAATAAATAACTTGTTAGAGAGCAATGTAAGGTTAAGGGAAAGCTTTCAGGTTTGTTTTGAAGAACGAGAAATACCAAATGGTGCTTGCAAGCAATGAGAAA".getBytes()));
            Assert.assertEquals(alignments, Arrays.asList(
                    Collections.singletonList(new MiniMap2Alignment(0, 0, false, 60, "120M")),
                    Collections.singletonList(new MiniMap2Alignment(0, 8, false, 60, "8S112M")),
                    Collections.singletonList(new MiniMap2Alignment(0, 0, false, 60, "112M8S")),
                    Collections.singletonList(new MiniMap2Alignment(0, 8, false, 60, "8S104M8S")),
                    Collections.singletonList(new MiniMap2Alignment(1, 60, true, 60, "120M"))));
        }
    }
}
