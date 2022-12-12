package dat.colourschemes;

import java.util.HashMap;
import java.util.Map;


final public class Hydro{


    private static final Map<String, String> hydro;
    static {
        hydro = new HashMap<String, String>();
        hydro.put("A", "#AD0052");
        hydro.put("R", "#0000FF");
        hydro.put("N", "#0C00F3");
        hydro.put("D", "#0C00F3");
        hydro.put("C", "#C2003D");
        hydro.put("Q", "#0C00F3");
        hydro.put("E", "#0C00F3");
        hydro.put("G", "#6A0095");
        hydro.put("H", "#1500EA");
        hydro.put("I", "#FF0000");
        hydro.put("L", "#EA0015");
        hydro.put("K", "#0000FF");
        hydro.put("M", "#B0004F");
        hydro.put("F", "#CB0034");
        hydro.put("P", "#4600B9");
        hydro.put("S", "#5E00A1");
        hydro.put("T", "#61009E");
        hydro.put("W", "#5B00A4");
        hydro.put("Y", "#4F00B0");
        hydro.put("V", "#F60009");
        hydro.put("B", "#0C00F3");
        hydro.put("X", "#680097");
        hydro.put("Z", "#0C00F3");
    }



    public static String getColour(char base) {
        return hydro.get(Character.toString(base));
    }
}