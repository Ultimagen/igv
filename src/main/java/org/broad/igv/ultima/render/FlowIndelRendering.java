package org.broad.igv.ultima.render;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.sam.*;
import org.broad.igv.track.RenderContext;

import java.awt.*;

public class FlowIndelRendering {

    // constants
    private static final String TAG_T0 = "t0";
    private static final String ATTR_TP = "tp";
    private static final String RG_ATTR_PL = "PL";
    private static final String RG_ATTR_MC = "mc";
    private static final String RG_ATTR_PL_ULTIMA = "ULTIMA";
    private static final double MIN_POSSIBLE_QUALITY = 0;

    // the color map
    private static ColorMap indelColorMap = ColorMap.getJet(42);

    // an Hmer
    static class Hmer {
        int start;
        int end;

        int size() { return end - start + 1; }
    }


    // is this alignment handled by this renderer?
    public boolean handlesAlignment(final Alignment alignment) {

        // we only handle sam alignments
        if ( !(alignment instanceof SAMAlignment) )
            return false;
        final SAMAlignment samAlignment = (SAMAlignment)alignment;

        // must be a modern Ultima flow
        final SAMReadGroupRecord readGroup = samAlignment.getRecord().getReadGroup();
        if ( !RG_ATTR_PL_ULTIMA.equals(readGroup.getAttribute(RG_ATTR_PL))
             &&  (readGroup.getAttribute(RG_ATTR_MC) == null) )
            return false;
        if ( !samAlignment.getRecord().hasAttribute(ATTR_TP)  )
            return false;

        return true;
    }

    // is this block (INS) handled (accepted)
    public boolean handlesBlock(AlignmentBlock block) {
        return true;
    }

    // is this gap (DEL) handled (accepted)
    public boolean handlesGap(Gap gap) {
        return true;
    }

    public void renderSmallInsertion(Alignment alignment,
                                     AlignmentBlock block,
                                     RenderContext context,
                                     int h, int x, int y,
                                     AlignmentTrack.RenderOptions renderOptions) {

        int pxWing = (h > 10 ? 2 : (h > 5) ? 1 : 0);
        int hairline = 2;
        if ( renderOptions.isInsertQualColoring() ) {
            pxWing = Math.min(pxWing, Math.max(1, (int) (1 / context.getScale())));
            hairline = Math.min(hairline, pxWing);
        }
        Graphics2D g = context.getGraphics();
        g.fillRect(x, y, hairline, h);
        g.fillRect(x - pxWing, y, hairline + 2 * pxWing, hairline);
        g.fillRect(x - pxWing, y + h - hairline, hairline + 2 * pxWing, hairline);


        // draw
        double q = getInsertionQuality(alignment, block);
        if ( !Double.isNaN(q) ) {
            Color currentColor = g.getColor();
            g.setColor(new Color(indelColorMap.getColor((int) q)));
            g.fillRect(x - pxWing, (int) (y + (h - hairline) * (q / 42)) - 1, hairline + 2 * pxWing, hairline * 2);
            g.setColor(currentColor);
        }
    }

    public void renderSmallInsertionWings(Alignment alignment,
                                     AlignmentBlock block,
                                     RenderContext context,
                                     int pxH, int pxTop, int pxRight, int pxLeft,
                                     AlignmentTrack.RenderOptions renderOptions) {

        int pxWing = (pxH > 10 ? 2 : 1);  // width of the cursor "wing"
        Graphics2D g = context.getGraphics();

        // adjust wing and hairline
        int hairline = 2;
        double locScale = context.getScale();
        if ( renderOptions.isInsertQualColoring() ) {
            pxWing = Math.min(pxWing, Math.max(1, (int) (1 / locScale)));
            hairline = Math.min(hairline, pxWing);
        }

        Color currentColor = g.getColor();
        g.setColor(AlignmentRenderer.purple);
        g.fillRect(pxLeft - pxWing, pxTop, pxRight - pxLeft + hairline * pxWing, hairline);
        g.fillRect(pxLeft - pxWing, pxTop + pxH - hairline, pxRight - pxLeft + hairline * pxWing, hairline);
        g.setColor(currentColor);

        // draw
        double q = getInsertionQuality(alignment, block);
        if ( !Double.isNaN(q) ) {
            currentColor = g.getColor();
            g.setColor(new Color(indelColorMap.getColor((int) q)));
            g.fillRect(pxLeft - pxWing, (int) (pxTop + (pxH - hairline) * (q / 42)), pxRight - pxLeft + hairline * pxWing, hairline);
            g.setColor(currentColor);
        }
    }

    public void renderDeletionGap(Alignment alignment,
                                          Gap gap,
                                          int y, int h, int x, int w,
                                          RenderContext context,
                                          AlignmentTrack.RenderOptions renderOptions) {

        // establish alignment blocks wrapping this gap (before and after)
        AlignmentBlock[] blocks = getGapWrappingBlocks(alignment, gap);
        if ( blocks == null )
            return;
        AlignmentBlock abPrev = blocks[0];
        AlignmentBlock abNext = blocks[1];

        // establish the nature of the gap
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        char        nextBlcokFirstBase = Character.toUpperCase((char)abNext.getBases().getByte(0));
        char        gapLastBase = Character.toUpperCase((char)genome.getReference(alignment.getChr(), gap.getStart() + gap.getnBases() - 1));
        boolean     isHmer = gapLastBase == nextBlcokFirstBase;

        // figure out quality to plot - if any
        SAMRecord record = ((SAMAlignment)alignment).getRecord();
        double q = Double.NaN;
        if ( isHmer ) {
            Hmer hmer = findHmer(record, abNext.getIndexOnRead(), (byte)gapLastBase, true, true);
            if ( hmer.size() >= getMC(record) ) {
                // HMER - length is at least max-hmer
                q = MIN_POSSIBLE_QUALITY;
            } else {
                // HMER - otherwise try TP
                q = getQualityFromTP(record, hmer, gap.getnBases());
            }
        } else {
            if ( gap.getnBases() == 1 ) {
                // NON-HMER, single base, try T0
                double qBefore = getQualityFromT0(record, abPrev, false);
                double qAfter = getQualityFromT0(record, abNext, true);
                if ( Double.isNaN(qBefore) )
                    q = qAfter;
                else if ( Double.isNaN(qAfter) )
                    q = qBefore;
                else
                    q = Math.max(qBefore, qAfter);
            }
        }
        if ( Double.isNaN(q) )
            return;

        // draw the marker
        int hairline = Math.min(2, (int) (1 / context.getScale()));
        Graphics2D g = context.getGraphics();
        Color currentColor = g.getColor();
        g.setColor(new Color(indelColorMap.getColor((int)q)));
        g.fillRect(x + (int) (w * q / 42) - hairline, y, hairline * 2, h);
        g.setColor(currentColor);
    }

    private double getInsertionQuality(Alignment alignment, AlignmentBlock block) {

        // establish alignment blocks wrapping this insertion block (before and after)
        AlignmentBlock[] blocks = getBlockWrappingBlocks(alignment, block);
        if ( blocks == null )
            return Double.NaN;
        AlignmentBlock abPrev = blocks[0];
        AlignmentBlock abNext = blocks[1];

        // establish the nature of the block
        char        nextBlcokFirstBase = Character.toUpperCase((char)abNext.getBases().getByte(0));
        char        blockLastBase = Character.toUpperCase((char)block.getBases().getByte(block.getBasesLength() - 1));
        boolean     isForwardHmer = blockLastBase == nextBlcokFirstBase;
        char        prevBlcokLastBase = Character.toUpperCase((char)abPrev.getBases().getByte(abPrev.getBasesLength() - 1));
        char        blockFirstBase = Character.toUpperCase((char)block.getBases().getByte(0));
        boolean     isBackwardsHmer = blockFirstBase == prevBlcokLastBase;

        // if not an hmer - nothing to do here
        if ( !isBackwardsHmer && !isForwardHmer )
            return Double.NaN;

        SAMRecord record = ((SAMAlignment)alignment).getRecord();
        Double forwardQ = Double.NaN;
        if ( isForwardHmer ) {
            Hmer hmer = findHmer(record, abNext.getIndexOnRead(), (byte)nextBlcokFirstBase, true, true);
            Hmer insertionHmer = findHmer(record, abNext.getIndexOnRead() - 1, (byte)nextBlcokFirstBase, true, false);
            forwardQ = getQualityFromTP(record, hmer, -insertionHmer.size());
        }
        Double backwardsQ = Double.NaN;
        if ( isBackwardsHmer ) {
            Hmer hmer = findHmer(record, abPrev.getIndexOnRead() + abPrev.getBasesLength() - 1 , (byte)prevBlcokLastBase, true, true);
            Hmer insertionHmer = findHmer(record, abPrev.getIndexOnRead() + abPrev.getBasesLength(), (byte)prevBlcokLastBase, false, true);
            backwardsQ = getQualityFromTP(record, hmer, -insertionHmer.size());
        }

        // integrate values
        if ( Double.isNaN(forwardQ) )
            return backwardsQ;
        else if ( Double.isNaN(backwardsQ) )
            return forwardQ;
        else
            return Math.max(backwardsQ, forwardQ);
    }

    private AlignmentBlock[] getBlockWrappingBlocks(Alignment alignment, AlignmentBlock block) {

        AlignmentBlock abPrev = null;
        AlignmentBlock abNext = null;
        for ( final AlignmentBlock b : alignment.getAlignmentBlocks() ) {
            if ( b.getEnd()  == block.getStart() )
                abPrev = b;
            else if ( b.getStart() == block.getStart() )
                abNext = b;
        }

        if ( abPrev != null && abNext != null ) {
            return new AlignmentBlock[] {abPrev, abNext};
        } else {
            return null;
        }
    }

    private AlignmentBlock[] getGapWrappingBlocks(Alignment alignment, Gap gap) {

        boolean blockFound = false;
        int blockIndex = 0;
        final AlignmentBlock[] blocks = alignment.getAlignmentBlocks();
        for ( ; blockIndex < blocks.length ; blockIndex++ ) {
            if ( blocks[blockIndex].getEnd() == gap.getStart() ) {
                blockFound = true;
                break;
            }
            else if ( blocks[blockIndex].getEnd() > gap.getStart() )
                break;
        }

        if ( blockFound && (blockIndex < blocks.length - 1) ) {
            return new AlignmentBlock[] {alignment.getAlignmentBlocks()[blockIndex], alignment.getAlignmentBlocks()[blockIndex + 1]};
        } else {
            return null;
        }
    }

    private double getQualityFromTP(SAMRecord record, Hmer hmer, int tpValue) {

        // get quals and tp
        byte[]      tp = record.getByteArrayAttribute(ATTR_TP);

        // scan for tpValue,
        for ( int ofs = hmer.start ; ofs <= hmer.end ; ofs++ ) {
            if ( tp[ofs] == tpValue ) {
                return record.getBaseQualities()[ofs];
            }
        }

        // if here, not found
        return Double.NaN;
    }

    private double getQualityFromT0(SAMRecord record, AlignmentBlock block, boolean delIsBeforeBlock) {

        final byte[] t0 = record.getStringAttribute(TAG_T0).getBytes();
        final int t0Index = delIsBeforeBlock ? block.getIndexOnRead() : (block.getIndexOnRead() + block.getLength() - 1);
        return t0[t0Index] - 33;
    }

    private int getMC(SAMRecord record) {
        try {
            return Integer.parseInt(record.getReadGroup().getAttribute("mc"));
        } catch (Exception e) {
            return 0;
        }
    }

    private Hmer findHmer(SAMRecord record, int start, byte base, boolean walkBackwards, boolean walkForward) {

        Hmer hmer = new Hmer();
        hmer.end = hmer.start = start;
        byte[] bases = record.getReadBases();
        if ( walkBackwards ) {
            while (hmer.start > 0 && bases[hmer.start - 1] == base)
                hmer.start--;
        }
        if ( walkForward ) {
            while ((hmer.end + 1) < bases.length && bases[hmer.end + 1] == base)
                hmer.end++;
        }

        return hmer;
    }


}
