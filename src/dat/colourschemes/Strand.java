package dat.colourschemes;

import java.util.HashMap;
import java.util.Map;


final public class Strand{


    private static final Map<String, String> strand;
    static {
        strand = new HashMap<String, String>();
        strand.put("A", "#5858A7");
        strand.put("R", "#6B6B94");
        strand.put("N", "#64649B");
        strand.put("D", "#2121DE");
        strand.put("C", "#9D9D62");
        strand.put("Q", "#8C8C73");
        strand.put("E", "#0000FF");
        strand.put("G", "#4949B6");
        strand.put("H", "#60609F");
        strand.put("I", "#ECEC13");
        strand.put("L", "#B2B24D");
        strand.put("K", "#4747B8");
        strand.put("M", "#82827D");
        strand.put("F", "#C2C23D");
        strand.put("P", "#2323DC");
        strand.put("S", "#4949B6");
        strand.put("T", "#9D9D62");
        strand.put("W", "#C0C03F");
        strand.put("Y", "#D3D32C");
        strand.put("V", "#FFFF00");
        strand.put("B", "#4343BC");
        strand.put("X", "#797986");
        strand.put("Z", "#4747B8");
    }



    public static String getColour(char base) {
        return strand.get(Character.toString(base));
    }
}