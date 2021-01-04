package asr;

import dat.phylo.TreeDecor;
import dat.phylo.TreeInstance;

import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class ThreadedDecorators<E> {
    final private ExecutorService executor;
    private Map<TreeDecor<E>, ThreadedDecorator<E>> jobs;

    public ThreadedDecorators(int nThreads) {
        executor = Executors.newFixedThreadPool(nThreads);
        jobs = new HashMap<>();
    }
    public TreeDecor<E> addDecorator(TreeDecor<E> decorator) {
        ThreadedDecorator<E> job = new ThreadedDecorator<>(decorator);
        jobs.put(decorator, job);
        return job;
    }

    public void runJobs() {
        for (ThreadedDecorator<E> job : jobs.values()) {
            if (!job.isCompleted())
                executor.execute(job);
        }
        executor.shutdown();
        while (!executor.isTerminated()) {}
    }

    public class ThreadedDecorator<E> implements Runnable, TreeDecor<E> {

        private TreeDecor<E> decorator;
        private boolean completed = false; // finished job?
        private TreeInstance ti;

        public ThreadedDecorator(TreeDecor<E> decorator) {
            this.decorator = decorator;
        }

        public E getDecoration(int idx) {
            if (!completed)
                runJobs();
            return decorator.getDecoration(idx);
        }

        public void decorate(TreeInstance ti) {
            this.ti = ti;
            this.completed = false;
        }
        public boolean isCompleted() {
            return completed;
        }
        @Override
        public void run() {
            decorator.decorate(ti);
            this.completed = true;
        }
    }

}
