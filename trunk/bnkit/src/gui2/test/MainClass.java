package gui2.test;


import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.awt.dnd.DnDConstants;
import java.awt.dnd.DragGestureEvent;
import java.awt.dnd.DragGestureListener;
import java.awt.dnd.DragSource;
import java.awt.dnd.DragSourceContext;
import java.awt.dnd.DragSourceDragEvent;
import java.awt.dnd.DragSourceDropEvent;
import java.awt.dnd.DragSourceEvent;
import java.awt.dnd.DragSourceListener;
import java.awt.dnd.DropTarget;
import java.awt.dnd.DropTargetDragEvent;
import java.awt.dnd.DropTargetDropEvent;
import java.awt.dnd.DropTargetEvent;
import java.awt.dnd.DropTargetListener;
import java.awt.dnd.InvalidDnDOperationException;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;

public class MainClass extends JFrame implements ActionListener, DropTargetListener {
  DragLabel source = new DragLabel("Drag and drop me to the following JButton", JLabel.CENTER);

  JButton target = new JButton();

  MainClass(String title) {
    super(title);
    setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    source.setForeground(Color.red);
    getContentPane().add(source, BorderLayout.NORTH);

    target.addActionListener(this);
    getContentPane().add(target, BorderLayout.SOUTH);

    new DropTarget(target, DnDConstants.ACTION_COPY_OR_MOVE, this);

    setSize(205, 100);

    setVisible(true);
  }

  public void actionPerformed(ActionEvent e) {
    JButton b = (JButton) e.getSource();
    b.setText("");
    source.setText("Drag and drop me to the following JButton");
  }

  public void dragEnter(DropTargetDragEvent e) {
    System.out.println("Entering drop target #1");
  }

  public void dragExit(DropTargetEvent e) {
    System.out.println("Exiting drop target #1");
  }

  public void dragOver(DropTargetDragEvent e) {
    System.out.println("Dragging over drop target #1");
  }

  public void drop(DropTargetDropEvent e) {
    System.out.println("Dropping");

    try {
      Transferable t = e.getTransferable();

      if (e.isDataFlavorSupported(DataFlavor.stringFlavor)) {
        e.acceptDrop(e.getDropAction());

        String s;
        s = (String) t.getTransferData(DataFlavor.stringFlavor);

        target.setText(s);

        e.dropComplete(true);
      } else
        e.rejectDrop();
    } catch (java.io.IOException e2) {
    } catch (UnsupportedFlavorException e2) {
    }
  }

  public void dropActionChanged(DropTargetDragEvent e) {
    System.out.println("Drop action changed #1");
  }

  public static void main(String[] args) {
    new MainClass("Drag and Drop Demo");
  }
}

class DragLabel extends JLabel implements DragGestureListener, DragSourceListener {
  private DragSource ds = DragSource.getDefaultDragSource();

  public DragLabel(String s, int alignment) {
    super(s, alignment);

    int action = DnDConstants.ACTION_COPY_OR_MOVE;
    ds.createDefaultDragGestureRecognizer(this, action, this);
  }

  public void dragGestureRecognized(DragGestureEvent e) {
    try {
      Transferable t = new StringSelection(getText());

      e.startDrag(DragSource.DefaultCopyNoDrop, t, this);
    } catch (InvalidDnDOperationException e2) {
      System.out.println(e2);
    }
  }

  public void dragDropEnd(DragSourceDropEvent e) {
    System.out.println("Drag and drop end");

    if (e.getDropSuccess() == false) {
      System.out.println("unsuccessful");
      return;
    }

    int action = e.getDropAction();
    if ((action & DnDConstants.ACTION_MOVE) != 0)
      setText("");
  }

  public void dragEnter(DragSourceDragEvent e) {
    System.out.println("Entering drop target #2");

    DragSourceContext ctx = e.getDragSourceContext();

    int action = e.getDropAction();
    if ((action & DnDConstants.ACTION_COPY) != 0)
      ctx.setCursor(DragSource.DefaultCopyDrop);
    else
      ctx.setCursor(DragSource.DefaultCopyNoDrop);
  }

  public void dragExit(DragSourceEvent e) {
    System.out.println("Exiting drop target #2");
  }

  public void dragOver(DragSourceDragEvent e) {
    System.out.println("Dragging over drop target #2");
  }

  public void dropActionChanged(DragSourceDragEvent e) {
    System.out.println("Drop action changed #2");
  }
}

           
       