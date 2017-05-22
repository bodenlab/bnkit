package dat.colourschemes;

import java.util.HashMap;
import java.util.Map;


final public class Clustal {


    private static final Map<String, String> clustal;
        static {
            clustal = new HashMap<String, String>();
            clustal.put("A", "#80B3E6");
            clustal.put("R", "#E6331A");
            clustal.put("N", "#1ACC1A");
            clustal.put("D", "#CC4DCC");
            clustal.put("C", "#E68080");
            clustal.put("Q", "#1ACC1A");
            clustal.put("E", "#CC4DCC");
            clustal.put("G", "#E6994D");
            clustal.put("H", "#1AB3B3");
            clustal.put("I", "#80B3E6");
            clustal.put("L", "#80B3E6");
            clustal.put("K", "#E6331A");
            clustal.put("M", "#80B3E6");
            clustal.put("F", "#80B3E6");
            clustal.put("P", "#CCCC00");
            clustal.put("S", "#1ACC1A");
            clustal.put("T", "#1ACC1A");
            clustal.put("W", "#80B3E6");
            clustal.put("Y", "#1AB3B3");
            clustal.put("V", "#80B3E6");
    }



    public static String getColour(char base) {

        return clustal.get(Character.toString(base));
    }
}