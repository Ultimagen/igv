package org.broad.igv.ultima.render;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.sam.*;
import org.broad.igv.track.RenderContext;

import java.awt.*;
import java.util.concurrent.atomic.AtomicBoolean;

public class FlowIndelRendering {

    public static final String TAG_T0 = "t0";
    private static final boolean GARRETY_LOW_PROB_MODE = true;
    private static ColorMap indelColorMap = ColorMap.getJet(42);
    private static final double MIN_PROB_DEFAULT = 0.01;
    private static final double MAX_PROB_DEFAULT = 0.9999;

    private static String ATTR_TP = "tp";
    private static String ATTR_TI = "ti";
    private static String RG_ATTR_PL = "PL";

    private static String RG_ATTR_MC = "mc";
    private static String RG_ATTR_PL_ULTIMA = "ULTIMA";

    public boolean handlesAlignment(final Alignment alignment) {

        if ( !(alignment instanceof SAMAlignment) )
            return false;
        final SAMAlignment samAlignment = (SAMAlignment)alignment;

        // must be a flow
        return getUltimaFileVersion(alignment) == UltimaFileFormat.BASE_TP;
    }

    public boolean handlesBlock(AlignmentBlock block) {
        return block.getLength() < 2;
    }

    public boolean handlesGap(Gap gap) {
        return true; //return gap.getnBases() < 2;
    }

    public void renderSmallInsertion(Alignment alignment,
                                     AlignmentBlock aBlock,
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

        double p = qualsAsProb(aBlock.getQualities());
        if ( p != 0 ) {
            double q = -10 * Math.log10(p);
            Color qColor = new Color(indelColorMap.getColor((int) q));
            Color c = g.getColor();
            g.setColor(qColor);
            g.fillRect(x - pxWing, (int) (y + (h - hairline) * (q / 42)) - 1, hairline + 2 * pxWing, hairline * 2);
            g.setColor(c);
        }
    }

    public void renderSmallInsertionWings(Alignment alignment,
                                     AlignmentBlock insertionBlock,
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

        Color c = g.getColor();
        g.setColor(AlignmentRenderer.purple);
        g.fillRect(pxLeft - pxWing, pxTop, pxRight - pxLeft + hairline * pxWing, hairline);
        g.fillRect(pxLeft - pxWing, pxTop + pxH - hairline, pxRight - pxLeft + hairline * pxWing, hairline);
        g.setColor(c);

        // Ultima: large (>1) INSERT case
        // map qual into a sort of a linear scale and add indicator
        double p = qualsAsProb(insertionBlock.getQualities());
        if ( p != 0 ) {
            double q = -10 * Math.log10(p);
            Color qColor = new Color(indelColorMap.getColor((int) q));
            c = g.getColor();
            g.setColor(qColor);
            g.fillRect(pxLeft - pxWing, (int) (pxTop + (pxH - hairline) * (q / 42)), pxRight - pxLeft + hairline * pxWing, hairline);
            g.setColor(c);
        }
    }

    public void renderDeletionGap(Alignment alignment,
                                          Gap gap,
                                          int y, int h, int x, int w,
                                          RenderContext context,
                                          AlignmentTrack.RenderOptions renderOptions) {

        // collect quals (experimental)
        Color[]       markerColor = new Color[2];
        double[]      markerQ = new double[2];
        AtomicBoolean markerFromT0_0 = new AtomicBoolean();
        AtomicBoolean markerFromT0_1 = new AtomicBoolean();

        // locate block who's end is the same as the start of the gap
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
            AlignmentBlock abPrev = alignment.getAlignmentBlocks()[blockIndex];
            AlignmentBlock abNext = alignment.getAlignmentBlocks()[blockIndex + 1];
            if (abPrev.getQualities().length > 0 && abNext.getQualities().length > 0) {

                // calc based on reference
                Genome genome = GenomeManager.getInstance().getCurrentGenome();
                char        gapBase0 = Character.toUpperCase((char)genome.getReference(alignment.getChr(), gap.getStart()));
                char        gapBase1 = Character.toUpperCase((char)genome.getReference(alignment.getChr(), gap.getStart() + gap.getnBases() - 1));
                char        alignBase0 = Character.toUpperCase((char)abPrev.getBases().getByte(abPrev.getBases().length - 1));
                char        alignBase1 = Character.toUpperCase((char)abNext.getBases().getByte(0));
                char        gapBase0p = Character.toUpperCase((char)genome.getReference(alignment.getChr(), gap.getStart() - 1));
                char        gapBase1n = Character.toUpperCase((char)genome.getReference(alignment.getChr(), gap.getStart() + gap.getnBases()));
                boolean     snp0 = gapBase0p != alignBase0;
                boolean     snp1 = gapBase1n != alignBase1;

                if ( !snp0 ) {
                    int         delLength = 1;
                    while ( (delLength + 1) <= gap.getnBases() &&
                            (gapBase0 == Character.toUpperCase(genome.getReference(alignment.getChr(), gap.getStart() + delLength))) )
                        delLength++;

                    double p = qualsAsProbDeleteTP(((SAMAlignment) alignment), abPrev, delLength, false, gap, gapBase0p == gapBase0, markerFromT0_0);

                    if ( p != 0 ) {
                        markerQ[0] = -10 * Math.log10(p);
                        markerColor[0] = new Color(indelColorMap.getColor((int) markerQ[0]));
                    }
                }
                if ( !snp1 ) {
                    int         delLength = 1;
                    while ( (delLength + 1) <= gap.getnBases() &&
                            (gapBase1 == Character.toUpperCase(genome.getReference(alignment.getChr(), gap.getStart() + gap.getnBases() - delLength))) )
                        delLength++;

                    double p = qualsAsProbDeleteTP((SAMAlignment) alignment, abNext, delLength, true, gap, gapBase1n == gapBase1, markerFromT0_1);

                    if ( p != 0 ) {
                        markerQ[1] = -10 * Math.log10(p);
                        markerColor[1] = new Color(indelColorMap.getColor((int) markerQ[1]));
                    }
                }
            }

            // perform T0 priority
            if ( markerFromT0_0.get() && markerFromT0_1.get() ) {

                // if both, retain higher quality in slot 0 and erase slot 1.
                // if just 1 then move to 0
                if ( markerQ[1] > markerQ[0] ) {
                    markerQ[0] = markerQ[1];    // copy slot 1 into slot 0
                    markerColor[0] = markerColor[1];
                }
                markerColor[1] = null; // mark slot 1 as not used
            } else {
                // T0 is considered only if on both sides. if not, cancel it out
                if ( markerFromT0_0.get() )
                    markerColor[0] = null;
                if ( markerFromT0_1.get() )
                    markerColor[1] = null;
            }

            // draw delete markers
            Graphics2D g = context.getGraphics();
            if ( markerColor[0] != null || markerColor[1] != null ) {

                int hairline = Math.min(2, (int) (1 / context.getScale()));

                Color c = g.getColor();
                if ((markerQ[0] == markerQ[1]) || (gap.getnBases() == 1) || (markerColor[0] == null ^ markerColor[1] == null)) {

                    // draw a full line, average as needed
                    double q1 = markerQ[0];
                    Color c1 = markerColor[0];
                    if (c1 == null) {
                        q1 = markerQ[1];
                        c1 = markerColor[1];
                    } else if (markerColor[1] != null && markerQ[1] != q1) {
                        q1 = (q1 + markerQ[1]) / 2;
                        c1 = new Color(indelColorMap.getColor((int) q1));
                    }
                    g.setColor(c1);
                    g.fillRect(x + (int) (w * q1 / 42) - hairline, y, hairline * 2, h);
                } else {
                    if (markerColor[0] != null) {
                        g.setColor(markerColor[0]);
                        g.fillRect(x + (int) (w * markerQ[0] / 42) - hairline, y, hairline * 2, h * 3 / 4);
                    }
                    if (markerColor[1] != null) {
                        g.setColor(markerColor[1]);
                        g.fillRect(x + (int) (w * markerQ[1] / 42) - hairline, y + h / 4, hairline * 2, h * 3 / 4);
                    }
                }

                g.setColor(c);
            }
        }
    }

   private enum UltimaFileFormat {
        NON_FLOW,
        BASE_TI,
        BASE_TP
    };

    static private UltimaFileFormat getUltimaFileVersion(Alignment alignment)
    {
        // check for preconditions for a flow file.  note that we extend only SAMAlignment instances
        if ( !(alignment instanceof SAMAlignment) ) {
            return UltimaFileFormat.NON_FLOW;
        }
        final SAMAlignment samAlignment = (SAMAlignment)alignment;
        if ( !isUltimaFlowReadGroup(samAlignment.getRecord().getReadGroup()) ) {
            return UltimaFileFormat.NON_FLOW;
        }

        if ( alignment.getAttribute(ATTR_TI) != null )
            return UltimaFileFormat.BASE_TI;
        else if ( alignment.getAttribute(ATTR_TP) != null )
            return UltimaFileFormat.BASE_TP;
        else
            return UltimaFileFormat.NON_FLOW;
    }

    private static boolean isUltimaFlowReadGroup(SAMReadGroupRecord readGroup) {

        // must have a read group to begin with
        if ( readGroup == null )
            return false;

        // modern files have an ultima platform
        if ( RG_ATTR_PL_ULTIMA.equals(readGroup.getAttribute(RG_ATTR_PL)) )
            return true;

        // fall back on the presence of an mc (max-class)
        if ( readGroup.getAttribute(RG_ATTR_MC) != null )
            return true;

        // if here, probably not an ultima flow read group
        return false;
    }

    static private double qualsAsProb(ByteSubarray quals) {

        // calc prob
        double              probSum = 0.0;
        double              probCount = 0;
        for ( int i = 0 ; i < quals.length ; i++ ) {
            final byte q = quals.getByte(i);
            if (q != 255) {
                probSum += Math.pow(10.0, -q / 10.0);
                probCount++;

            }
        }
        if ( probCount != 0 )
            return probSum / probCount;
        else
            return -1;
    }


    public double qualsAsProbInsertTP(SAMAlignment samAlignment, AlignmentBlock block) {
        return qualsAsProbInsertTP(samAlignment, block, 0, block.getLength());
    }

    private class HMer {
        int     start;
        int     end;
    }

    public double qualsAsProbInsertTP(SAMAlignment samAlignment, AlignmentBlock block, int fragOfs, int fragLength) {

        // for now, handle only hmer inserts
        if ( !blockIsHmer(block, fragOfs, fragLength) ) {

            // by definition, blocks failing this test can only be of multiple bases
            // try breaking them into two block fragments, one from front and one from back
            byte[]  bases = block.getBases().getBytes();
            int     f1ofs = 0;
            int     f1Length = 0;
            int     base = bases[f1ofs];
            for ( int n = 0 ; n < bases.length ; n++ ) {
                if (bases[n] == base)
                    f1Length++;
                else
                    break;
            }

            int     f2ofs = bases.length;
            int     f2Length = 0;
            base = bases[f2ofs-1];
            for ( int n = bases.length - 1 ; n >= 0 ; n-- ) {
                if (bases[n] == base) {
                    f2Length++;
                    f2ofs--;
                }
                else
                    break;
            }

            double      p1 = qualsAsProbInsertTP(samAlignment, block, f1ofs, f1Length);
            double      p2 = qualsAsProbInsertTP(samAlignment, block, f2ofs, f2Length);

            if ( p1 == 0 )
                return p2;
            else if ( p2 == 0 )
                return p1;
            else
                return (p1 + p2) / 2;
        }

        // access read/record
        SAMRecord record = samAlignment.getRecord();
        if ( record == null )
            return 0;


        // establish the hmer on which this block sits (inside the read...).
        byte        base = block.getBase(fragOfs);
        int         start = block.getIndexOnRead() + fragOfs;

        // short-circuit a simple and common case: insert of 1 with tp=-1 -> quality is already here!
        if ( fragLength == 1 ) {
            byte[]      tp = record.getByteArrayAttribute(ATTR_TP);
            if ( tp[start] == -1 ) {
                byte      q = block.getQuality(0);
                return Math.pow(10.0, -q / 10.0);
            }
        }
        HMer        hmer = findHmer(record, start, fragLength, base);

        // find tp value and return it
        return findQualByTPValue(record, hmer, -fragLength);
    }

    private HMer findHmer(SAMRecord record, int start, int length, byte base) {

        HMer        hmer = new HMer();
        hmer.start = start;

        byte[]      bases = record.getReadBases();
        hmer.end = hmer.start + (length - 1);
        while ( hmer.start > 0 && bases[hmer.start - 1] == base )
            hmer.start--;
        while ( (hmer.end + 1) < bases.length  && bases[hmer.end + 1] == base )
            hmer.end++;

        return hmer;
    }

    private double findQualByTPValue(SAMRecord record, HMer hmer, int tpValue) {

        // get quals and tp
        byte[]      tp = record.getByteArrayAttribute(ATTR_TP);

        // scan for tpValue, extract qual, choose min if more than one
        int foundQual = Integer.MAX_VALUE;
        boolean exactMatchFound = false;
        for ( int ofs = hmer.start ; ofs <= hmer.end ; ofs++ ) {
            boolean tpApplicable = false;
            if ( tp[ofs] == tpValue ) {
                exactMatchFound = true;
                tpApplicable = true;
            }
            if (GARRETY_LOW_PROB_MODE)
                if ( tpValue > 0 && tp[ofs] < tpValue )
                    tpApplicable = true;
            if ( tpApplicable ) {
                foundQual = Math.min(foundQual, record.getBaseQualities()[ofs]);
            }
        }

        // if here, none of the tp values matched exactly
        // check for the special case of a delete which is of an original hmer larger than mc (def:12)
        if ( !exactMatchFound && tpValue > 0 ) {
            int     mc = getMC(record);
            int     hmerSize = hmer.end - hmer.start + 1;
            if ( tpValue + hmerSize > mc )
                return 1.0 - MIN_PROB_DEFAULT;
        }

        // if we managed to get a qual value from other TP values, use them
        if ( foundQual != Integer.MAX_VALUE )
            return Math.pow(10.0, -foundQual / 10.0);

        // if here and a deletion, return minimal probability as well.
        // this is overlapping with the previous case, but as part of the visual development
        // it is considered a separate case
        if ( tpValue > 0)
            return 1 - MAX_PROB_DEFAULT;

        // if here, simply fail
        return 0;
    }

    private int getMC(SAMRecord record) {
        try {
            return Integer.parseInt(record.getReadGroup().getAttribute("mc"));
        } catch (Exception e) {
            return 0;
        }
    }

    private boolean blockIsHmer(AlignmentBlock block, int fragOfs, int fragLength) {
        byte[]      bases = block.getBases().getBytes();

        if ( bases == null || bases.length < (fragOfs + fragLength) )
            return false;
        else if ( fragLength == 1 )
            return true;

        byte base = bases[fragOfs];
        for ( int n = fragOfs + 1 ; n < fragOfs + fragLength ; n++ ) {
            if (bases[n] != base)
                return false;
        }
        return true;
    }

    public double qualsAsProbDeleteTP(SAMAlignment samAlignment, AlignmentBlock block, int delLength, boolean delIsBeforeBlock, Gap gap, boolean belongsToHmer, AtomicBoolean fromT0)
    {
        // access read/record
        SAMRecord   record = samAlignment.getRecord();
        if ( record == null )
            return 0;

        // establish the hmer on which this block sits (inside the read...)
        byte        base = block.getBase(delIsBeforeBlock ? 0 : block.getLength() - 1);
        int         start = block.getIndexOnRead();
        if ( !delIsBeforeBlock )
            start += (block.getLength() - 1);
        HMer        hmer = findHmer(record, start, 0, base);

        // try establishing by using t0
        final double t0qual = qualsAsProbDeleteTPByT0(samAlignment, record, block, delLength, delIsBeforeBlock, gap);
        if ( t0qual != 0 ) {
            fromT0.set(true);
            return t0qual;
        }

        return belongsToHmer ? findQualByTPValue(record, hmer, delLength) : 0;
    }

    private double qualsAsProbDeleteTPByT0(SAMAlignment samAlignment, SAMRecord record, AlignmentBlock block, int delLength, boolean delIsBeforeBlock, Gap gap) {

        // consider using t0 only if DEL of one base (additional conditions to follow)
        if (delLength != 1) {
            return 0;
        }

        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        final byte delBase = (byte) Character.toUpperCase(genome.getReference(samAlignment.getChr(), gap.getStart()));

        // delete is before the block?
        if (delIsBeforeBlock) {

            if (delBase == record.getReadBases()[block.getIndexOnRead()])
                return 0;

            // extract t0 value
            if (!record.hasAttribute(TAG_T0))
                return 0;
            final byte[] t0 = record.getStringAttribute(TAG_T0).getBytes();
            return Math.pow(10.0, (t0[block.getIndexOnRead()] - '!') / -10.0);

        } else {

            // delete is after the block
            if (delBase == record.getReadBases()[block.getIndexOnRead() + block.getLength() - 1])
                return 0;

            // extract t0 value
            if (!record.hasAttribute(TAG_T0))
                return 0;
            final byte[] t0 = record.getStringAttribute(TAG_T0).getBytes();
            return Math.pow(10.0, (t0[block.getIndexOnRead() + block.getLength() - 1] - '!') / -10.0);
        }
    }
}
