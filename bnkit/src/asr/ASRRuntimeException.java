package asr;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

public class ASRRuntimeException extends RuntimeException {

    private String msg = null;

    public ASRRuntimeException(String msg) {
        super(msg);
        this.msg = msg;
    }

    @Override
    public String getMessage() {
        if (msg != null)
            return msg;
        else
            return "Unknown ASR runtime exception";
    }
}
