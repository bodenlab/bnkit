/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package gui2;

/**
 *
 * @author Jun
 */

import java.util.Hashtable;

import javax.swing.JFrame;

import com.mxgraph.swing.mxGraphComponent;
import com.mxgraph.util.mxConstants;
import com.mxgraph.view.mxGraph;
import com.mxgraph.view.mxStylesheet;

public class StyleExample extends JFrame {

    private static final long serialVersionUID = 672772281200016954L;

    public static final String MY_CUSTOM_VERTEX_STYLE = "MY_CUSTOM_VERTEX_STYLE";
    public static final String MY_CUSTOM_EDGE_STYLE = "MY_CUSTOM_EDGE_STYLE";
    public static final String NEW_CUSTOM_EDGE = "NEW_CUSTOM_EDGE";

    private static void setStyleSheet(mxGraph graph) {

        Hashtable<String, Object> style;
        mxStylesheet stylesheet = graph.getStylesheet();

        // base style
        Hashtable<String, Object> baseStyle = new Hashtable<String, Object>();
        baseStyle.put(mxConstants.STYLE_STROKECOLOR, "#FF0000");

        // custom vertex style
        style = new Hashtable<String, Object>(baseStyle);
        style.put(mxConstants.STYLE_FILLCOLOR, "#FFFF00");
        stylesheet.putCellStyle(MY_CUSTOM_VERTEX_STYLE, style);

        // custom edge style
        style = new Hashtable<String, Object>(baseStyle);
        style.put(mxConstants.STYLE_STROKEWIDTH, 3);
        stylesheet.putCellStyle(MY_CUSTOM_EDGE_STYLE, style);

        // new custom vertex style -- Jun
        style = new Hashtable<String, Object>(baseStyle);
//        style.put(mxConstants.SHAPE_SWIMLANE, style);
        style.put(mxConstants.STYLE_SHAPE, mxConstants.SHAPE_ELLIPSE);
        stylesheet.putCellStyle(NEW_CUSTOM_EDGE, style);
        
        
        
    }

    public StyleExample() {
        super("Hello, World!");

        final mxGraph graph = new mxGraph();
        Object parent = graph.getDefaultParent();

        // create styles
        setStyleSheet(graph);

        graph.getModel().beginUpdate();
        try {

            Object v1 = graph.insertVertex(parent, null, "Hello", 20, 20, 80, 30, MY_CUSTOM_VERTEX_STYLE);
            Object v2 = graph.insertVertex(parent, null, "World!", 240, 150, 80, 30, MY_CUSTOM_VERTEX_STYLE);
            Object v3 = graph.insertVertex(parent, null, "text overfloooooooooww", 0, 0, 80, 30, NEW_CUSTOM_EDGE);
            graph.insertEdge(parent, null, "Edge", v1, v2, MY_CUSTOM_EDGE_STYLE);

        } finally {
            graph.getModel().endUpdate();
        }

        final mxGraphComponent graphComponent = new mxGraphComponent(graph);
        getContentPane().add(graphComponent);

    }

    public static void main(String[] args) {
        StyleExample frame = new StyleExample();
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(400, 320);
        frame.setVisible(true);
    }

}