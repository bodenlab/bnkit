package alignment.utilities;

/**
 * Int which can be incremented
 */
public class MutableInt {

    // Start at one because we are counting
    private int value = 1;

    public void increment() {
        ++value;
    }


    public int getValue() {
        return value;
    }

    public String toString(){
        return Integer.toString(value);
    }
}
