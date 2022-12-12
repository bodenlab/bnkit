package dat.colourschemes;

import java.util.HashMap;
import java.util.Map;


final public class Buried{


    private static final Map<String, String> buried;
    static {
        buried = new HashMap<String, String>();
        buried.put("A", "#00A35C");
        buried.put("R", "#00FC03");
        buried.put("N", "#00EB14");
        buried.put("D", "#00EB14");
        buried.put("C", "#0000FF");
        buried.put("Q", "#00F10E");
        buried.put("E", "#00F10E");
        buried.put("G", "#009D62");
        buried.put("H", "#00D52A");
        buried.put("I", "#0054AB");
        buried.put("L", "#007B84");
        buried.put("K", "#00FF00");
        buried.put("M", "#009768");
        buried.put("F", "#008778");
        buried.put("P", "#00E01F");
        buried.put("S", "#00D52A");
        buried.put("T", "#00DB24");
        buried.put("W", "#00A857");
        buried.put("Y", "#00E619");
        buried.put("V", "#005FA0");
        buried.put("B", "#00EB14");
        buried.put("X", "#00B649");
        buried.put("Z", "#00F10E");
    }



    public static String getColour(char base) {
        return buried.get(Character.toString(base));
    }
}