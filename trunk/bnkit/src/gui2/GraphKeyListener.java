/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gui2;

import com.mxgraph.model.mxIGraphModel;
import com.mxgraph.view.mxGraph;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;

/**
 * Handles key press events for graph
 *
 * @author jun
 */
public class GraphKeyListener extends KeyAdapter {

    mxGraph graph = null;

    mxIGraphModel model = null;

    public GraphKeyListener(mxGraph graph, mxIGraphModel model) {
        this.graph = graph;
        this.model = model;
    }

    public void keyPressed(KeyEvent event) {
        int key = event.getKeyCode();
        if (key == KeyEvent.VK_DELETE){
            System.out.println("Delete pressed");
        }
    }

}
