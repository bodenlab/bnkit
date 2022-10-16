package util;

import java.util.HashMap;
import java.util.Map;

public class MilliTimer {

    private Map<String, Long> tagmap = new HashMap<>();
    private Map<String, Long> elapsed = new HashMap<>();
    private Map<String, Integer> nstops = new HashMap<>();

    public void start(String tag) {
        tagmap.put(tag, System.currentTimeMillis());
    }

    public void stopStart(String stoptag, String starttag) {
        stop(stoptag, false);
        start(starttag);
    }
    public void stop(String tag) {
        stop(tag, false);
    }
    public void stop(String tag, boolean restart) {
        Long start = tagmap.get(tag);
        Long stop = System.currentTimeMillis();
        if (start != null) {
            Long sofar = elapsed.get(tag);
            Integer n = nstops.get(tag);
            if (sofar == null || n == null) {
                sofar = 0L;
                n = 0;
            }
            sofar += (stop - start);
            n += 1;
            elapsed.put(tag, sofar);
            nstops.put(tag, n);
        }
        if (restart)
            start(tag);
    }

    public void report() {
        report(false);
    }
    public void report(boolean total) {
        if (total)
            System.out.println("Tag    \tTotal ms");
        else
            System.out.println("Tag    \tMean ms");
        for (Map.Entry<String, Long> entry : elapsed.entrySet()) {
            Integer n = nstops.get(entry.getKey());
            System.out.println(entry.getKey() + "\t" + (total ? entry.getValue() : entry.getValue() / n));
        }
    }

}
