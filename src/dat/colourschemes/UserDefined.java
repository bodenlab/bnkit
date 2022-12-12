package dat.colourschemes;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;


final public class UserDefined {

//    FileReader fr = new FileReader();

    String dict = "K = red, N = green, P=blue";


    private static final Map<String, String> userDefined;
    static {
        userDefined = new HashMap<String, String>();
        userDefined.put("A", "white");
        userDefined.put("R", "white");
        userDefined.put("N", "white");
        userDefined.put("D", "white");
        userDefined.put("C", "white");
        userDefined.put("Q", "white");
        userDefined.put("E", "white");
        userDefined.put("G", "white");
        userDefined.put("H", "white");
        userDefined.put("I", "white");
        userDefined.put("L", "white");
        userDefined.put("K", "white");
        userDefined.put("M", "white");
        userDefined.put("F", "white");
        userDefined.put("P", "white");
        userDefined.put("S", "white");
        userDefined.put("T", "white");
        userDefined.put("W", "white");
        userDefined.put("Y", "white");
        userDefined.put("V", "white");
    }

//    public static void updateDict(String filename) throws FileNotFoundException, IOException{
//        FileReader fr = new FileReader(filename);
//        BufferedReader br = new BufferedReader(fr);
//
//        String line;
//
//        while((line = br.readLine()) != null){
//            String[] pairs = line.split(",");
//
//        }
//
//
//    }

    public static void updateDict(String filename) throws FileNotFoundException, IOException{
        String dict = "K = red, N = green, P=blue";

        String[] commaSplit = dict.split(",");







    }




    public static String getColour(char base) {

        return userDefined.get(Character.toString(base));
    }
}