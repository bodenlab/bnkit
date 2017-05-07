package dat.colourschemes;

import java.util.HashMap;
import java.util.Map;


final public class Turn{


    private static final Map<String, String> turn;
    static {
        turn = new HashMap<String, String>();
        turn.put("A", "#2CD3D3");
        turn.put("R", "#708F8F");
        turn.put("N", "#FF0000");
        turn.put("D", "#E81717");
        turn.put("C", "#A85757");
        turn.put("Q", "#3FC0C0");
        turn.put("E", "#778888");
        turn.put("G", "#FF0000");
        turn.put("H", "#708F8F");
        turn.put("I", "#00FFFF");
        turn.put("L", "#1CE3E3");
        turn.put("K", "#7E8181");
        turn.put("M", "#1EE1E1");
        turn.put("F", "#1EE1E1");
        turn.put("P", "#F60909");
        turn.put("S", "#E11E1E");
        turn.put("T", "#738C8C");
        turn.put("W", "#738C8C");
        turn.put("Y", "#9D6262");
        turn.put("V", "#07F8F8");
        turn.put("B", "#F30C0C");
        turn.put("X", "#7C8383");
        turn.put("Z", "#5BA4A4");
    }



    public static String getColour(char base) {
        return turn.get(Character.toString(base));
    }
}