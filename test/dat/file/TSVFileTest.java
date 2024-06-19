package dat.file;

import org.junit.jupiter.api.Test;

import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

class TSVFileTest {

    String[] examples = {
            "0.13;4-nitrophenyl alpha-D-maltoheptaoside-4,6-O-ethylidene_count=1",
            "377.0;starch_count=1;619.0;starch_count=2;880.0;starch_count=3",
            "7.2-7.5_count=1",
            "34_count=1;55_count=3",
            "SM00642;SM00632;",
            "None"};
    @Test
    void tokeniseByCount() {
        for (int i = 0; i < examples.length; i ++) {
            Object[] ys = TSVFile.tokeniseByCount(examples[i]);
            for (Object s : ys)
                System.out.print(s + "\t");
            System.out.println();
        }
    }
}