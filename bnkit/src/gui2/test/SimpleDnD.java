package gui2.test;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JTextField;
import javax.swing.TransferHandler;


public class SimpleDnD extends JFrame {

    JTextField field;
    JButton button;

    public SimpleDnD() {

        setTitle("Simple Drag & Drop");

        setLayout(null);

        button = new JButton("Button");
        button.setBounds(200, 50, 90, 25);

        field = new JTextField();
        field.setBounds(30, 50, 150, 25);

        add(button);
        add(field);

        field.setDragEnabled(true);
        button.setTransferHandler(new TransferHandler("text"));

        setSize(330, 150);
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setLocationRelativeTo(null);
        setVisible(true);
    }

    public static void main(String[] args) {
        new SimpleDnD();
    }
}