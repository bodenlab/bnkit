package bn.ctmc;
import bn.ctmc.matrix.JTTGap;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class GapSubstModelTest {


    @Test
    void testStationaryProbabilities() {
        SubstModel JTT  = SubstModel.createModel("JTT");
        GapSubstModel JTT_no_indel = new JTTGap(0, 0);

        for (Object state : JTT.getDomain().getValues()) {
            double probJTT = JTT.getProb(state);
            double probNoIndel = JTT_no_indel.getProb(state);
            assertEquals(probJTT, probNoIndel, 1e-10, "Stationary probabilities differ for state: " + state);
        }
    }

    @Test
    void testConditionalProbabilities() {
        SubstModel JTT = SubstModel.createModel("JTT");
        GapSubstModel JTT_no_indel = new JTTGap(0, 0);
        double t = 0.5;

        for (Object from : JTT.getDomain().getValues()) {
            for (Object to : JTT.getDomain().getValues()) {
                double probJTT = JTT.getProb(from, to, t);
                double probNoIndel = JTT_no_indel.getProb(from, to, t);
                assertEquals(probJTT, probNoIndel, 1e-5, String.format("Conditional probabilities differ for %s -> %s", from, to));
            }
        }
    }

    @Test
    void testGapProbabilities() {
        GapSubstModel JTT_no_indel = new JTTGap(0, 0);

        // Ensure that the gap probability is zero when there are no indels
        double gapProb = JTT_no_indel.getProb('-');
        assertEquals(0.0, gapProb, 1e-5, "Gap probability should be zero for JTT_no_indel");
    }

    @Test
    void testModelCopy() {
        GapSubstModel JTTGap = new JTTGap(0.05, 0.05);
        GapSubstModel JTTGapCopy = JTTGap.deepCopy();
        double t = 0.5;

        for (Object from : JTTGap.getDomain().getValues()) {
            for (Object to : JTTGap.getDomain().getValues()) {
                double probJTT = JTTGap.getProb(from, to, t);
                double probCopy = JTTGapCopy.getProb(from, to, t);
                assertEquals(probJTT, probCopy, 1e-5, String.format("Conditional probabilities differ for %s -> %s", from, to));
            }
        }
    }

}
