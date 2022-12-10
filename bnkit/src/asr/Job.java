package asr;

public interface Job {

    public String username = "anonymous";
    public String auth_token = null;

    public int NTHREADS = 4;

    public static
    int HOST_GENERIC = 1,    // e.g. single-user desktop, laptop etc
        HOST_WEBSERVER = 2,  // e.g. multi-user server with limited resources
        HOST_HPC = 3;        // e.g. multi-user server with plenty of resources

    public int HOST_ENVIRONMENT = HOST_GENERIC;

    public int PRIORITY = 0; // from -10 to +10, with -10 meaning lowest, +10 highest priority, 0 is neutral

}
