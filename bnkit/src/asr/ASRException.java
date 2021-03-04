package asr;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

public class ASRException extends Exception {

    private String msg = null;
    private FileIssues fileissues = null;

    public static class FileIssues {
        Map<String, Collection<String>> errors = new HashMap<>();
        public FileIssues() {
        }
        public FileIssues(String issue, Collection<String> params) {
            errors.put(issue, params);
        }
        public void add(String issue, Collection<String> params) {
            errors.put(issue, params);
        }
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            for (Map.Entry<String, Collection<String>> entry : errors.entrySet()) {
                sb.append(entry.getKey() + System.lineSeparator());
                int i = 0;
                for (String errvalue : entry.getValue())
                    sb.append(entry.getValue().toString() + (i ++ < entry.getValue().size() - 1 ? "; " : ""));
            }
            return sb.toString();
        }
        public boolean isLoaded() {
            return errors.size() > 0;
        }
    }

    public ASRException(String msg) {
        super(msg);
        this.msg = msg;
    }

    public ASRException(FileIssues issues) {
        super("ASR File issues");
        this.fileissues = issues;
    }

    @Override
    public String getMessage() {
        if (msg != null)
            return msg;
        else if (fileissues != null){
            return fileissues.toString();
        } else {
            return "ASR unknown exception";
        }
    }
}
