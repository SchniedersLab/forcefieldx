//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.ui;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.JViewport;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.ArrayList;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;

import ffx.utilities.Keyword;

/**
 * The KeywordComponent class is used to represent one TINKER keyword.
 *
 * @author Michael J. Schnieders
 */
public final class KeywordComponent implements MouseListener, ActionListener,
        ChangeListener, DocumentListener {

    private static final Logger logger = Logger.getLogger(KeywordComponent.class.getName());

    public enum SwingRepresentation {

        TEXTFIELD, CHECKBOX, CHECKBOXES, EDITCOMBOBOX, COMBOBOX, MULTIPOLE
    }

    /**
     * This is used to test if ANY keyword has been modified, so it is static
     * across all keyword objects
     */
    private static boolean isModified = false;
    private static Dimension labelDimension;
    private static Dimension entryDimension;

    static {
        JTextField textField = new JTextField(20);
        labelDimension = textField.getPreferredSize();
        textField = new JTextField(25);
        entryDimension = textField.getPreferredSize();
    }

    /**
     * TINKER Keyword String.
     */
    private String keyword;
    /**
     * TINKER Keyword Group.
     */
    private String keywordGroup;
    /**
     * An ArrayList of Components used to represent this Keyword
     */
    private ArrayList<Component> keywordValues;
    private JPanel keywordGUI = null;
    private final FlowLayout flowLayout;
    /**
     * The type of Swing Conponent used in representing this Keyword
     */
    private SwingRepresentation swingRepresentation;
    private String[] options;
    private String keywordDescription;
    private JTextArea output;
    private boolean active;
    private boolean init = false;

    /**
     * The Default Constructor k - Keyword String kg - Keyword Group t - Type of
     * GUI Components used to represent Keyword modifiers d - Keyword
     * description
     *
     * @param keyword             a {@link java.lang.String} object.
     * @param keywordGroup        a {@link java.lang.String} object.
     * @param swingRepresentation a {@link ffx.ui.KeywordComponent.SwingRepresentation} object.
     * @param keywordDescription  a {@link java.lang.String} object.
     * @param jTextArea           a {@link javax.swing.JTextArea} object.
     */
    KeywordComponent(String keyword, String keywordGroup, SwingRepresentation swingRepresentation,
                     String keywordDescription, JTextArea jTextArea) {
        flowLayout = new FlowLayout(FlowLayout.LEFT, 5, 5);
        this.keyword = keyword;
        this.keywordGroup = keywordGroup;
        keywordValues = new ArrayList<>();
        this.swingRepresentation = swingRepresentation;
        this.keywordDescription = keywordDescription;
        output = jTextArea;
        active = false;
        flowLayout.setHgap(2);
        flowLayout.setVgap(1);
    }

    /**
     * <p>
     * Constructor for KeywordComponent.</p>
     *
     * @param keyword             a {@link java.lang.String} object.
     * @param keywordGroup        a {@link java.lang.String} object.
     * @param swingRepresentation a {@link ffx.ui.KeywordComponent.SwingRepresentation} object.
     * @param keywordDescription  a {@link java.lang.String} object.
     * @param jTextArea           a {@link javax.swing.JTextArea} object.
     * @param options             an array of {@link java.lang.String} objects.
     */
    KeywordComponent(String keyword, String keywordGroup, SwingRepresentation swingRepresentation,
                     String keywordDescription, JTextArea jTextArea, String[] options) {
        this(keyword, keywordGroup, swingRepresentation, keywordDescription, jTextArea);
        this.options = options;
    }

    /**
     * <p>
     * fillPanel</p>
     *
     * @param jPanel             a {@link javax.swing.JPanel} object.
     * @param gridBagLayout      a {@link java.awt.GridBagLayout} object.
     * @param gridBagConstraints a {@link java.awt.GridBagConstraints} object.
     */
    static void fillPanel(JPanel jPanel, GridBagLayout gridBagLayout, GridBagConstraints gridBagConstraints) {
        JLabel jfill = new JLabel(" ");
        gridBagConstraints.weightx = 1;
        gridBagConstraints.weighty = 1;
        gridBagConstraints.fill = GridBagConstraints.BOTH;
        gridBagConstraints.gridwidth = GridBagConstraints.REMAINDER;
        gridBagConstraints.gridheight = GridBagConstraints.REMAINDER;
        gridBagLayout.setConstraints(jfill, gridBagConstraints);
        jPanel.add(jfill);
    }

    /**
     * <p>
     * isKeywordModified</p>
     *
     * @return a boolean.
     */
    static boolean isKeywordModified() {
        return isModified;
    }

    /**
     * <p>
     * setKeywordModified</p>
     *
     * @param b a boolean.
     */
    static void setKeywordModified(boolean b) {
        isModified = b;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    @SuppressWarnings("unchecked")
    public void actionPerformed(ActionEvent evt) {
        synchronized (this) {
            isModified = true;
            if (evt.getSource() instanceof JButton) {
                JButton button = (JButton) evt.getSource();
                if (button.getText().equalsIgnoreCase("Add")) {
                    JTextField text = (JTextField) keywordValues.get(3);
                    String s = text.getText();
                    if (s != null && !s.trim().equals("")) {
                        JComboBox<String> jcb = (JComboBox<String>) keywordValues.get(1);
                        jcb.addItem(s);
                        text.setText("");
                        active = true;
                    }
                } else if (button.getText().equalsIgnoreCase("Remove")) {
                    JComboBox<String> jcb = (JComboBox<String>) keywordValues.get(1);
                    if (jcb.getItemCount() > 0) {
                        jcb.removeItemAt(jcb.getSelectedIndex());
                    }
                    if (jcb.getItemCount() == 0) {
                        active = false;
                    }
                }
            } else if (evt.getSource() instanceof JComboBox) {
                JComboBox<String> jcb = (JComboBox<String>) evt.getSource();
                String selected = (String) jcb.getSelectedItem();
                if (selected == null) {
                    active = false;
                } else {
                    active = !selected.equalsIgnoreCase("DEFAULT");
                }
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void changedUpdate(DocumentEvent evt) {
        isModified = true;
    }

    private void checkBoxesInit() {
        checkBoxInit();
        for (String s : options) {
            checkBoxInit(s);
        }
    }

    private void checkBoxInit() {
        JCheckBox cb = new JCheckBox(keyword, false);
        cb.setPreferredSize(labelDimension);
        cb.setMaximumSize(labelDimension);
        cb.addMouseListener(this);
        cb.addChangeListener(this);
        keywordValues.add(cb);
    }

    private void checkBoxInit(String label) {
        JCheckBox cb = new JCheckBox(label, false);
        cb.addMouseListener(this);
        cb.addChangeListener(this);
        keywordValues.add(cb);
    }

    /**
     * <p>
     * clearKeywordComponent</p>
     */
    void clearKeywordComponent() {
        synchronized (this) {
            active = false;
            if (!init) {
                return;
            }
            if (swingRepresentation == SwingRepresentation.CHECKBOX) {
                ((JCheckBox) keywordValues.get(0)).setSelected(false);
            } else if (swingRepresentation == SwingRepresentation.CHECKBOXES) {
                for (Component keywordValue : keywordValues) {
                    ((JCheckBox) keywordValue).setSelected(false);
                }
            } else if (swingRepresentation == SwingRepresentation.COMBOBOX) {
                JComboBox jcb = (JComboBox) keywordValues.get(1);
                jcb.setSelectedItem("DEFAULT");
            } else if (swingRepresentation == SwingRepresentation.EDITCOMBOBOX) {
                ((JComboBox) keywordValues.get(1)).removeAllItems();
            } else if (swingRepresentation == SwingRepresentation.TEXTFIELD) {
                ((JTextField) keywordValues.get(1)).setText("");
            } else {
                logger.severe("Keyword Component: Unknown Keyword Type");
                logger.severe("Force Field X can not continue...");
                System.exit(-1);
            }
        }
    }

    private void comboBoxInit() {
        JLabel jl = new JLabel(keyword);
        jl.addMouseListener(this);
        jl.setPreferredSize(labelDimension);
        jl.setMaximumSize(labelDimension);
        keywordValues.add(jl);
        JComboBox<String> cb = new JComboBox<>();
        cb.setEditable(false);
        cb.addMouseListener(this);
        cb.addActionListener(this);
        cb.setPreferredSize(entryDimension);
        cb.setMaximumSize(entryDimension);
        for (String s : options) {
            cb.addItem(s);
        }
        cb.setSelectedItem("DEFAULT");
        keywordValues.add(cb);
    }

    private void editComboBoxInit() {
        JLabel jl = new JLabel(keyword);
        keywordValues.add(jl);
        jl.addMouseListener(this);
        jl.setPreferredSize(labelDimension);
        jl.setMaximumSize(labelDimension);
        JComboBox cb = new JComboBox();
        cb.setEditable(false);
        cb.addActionListener(this);
        cb.setPreferredSize(entryDimension);
        cb.setMaximumSize(entryDimension);
        keywordValues.add(cb);
        JButton remove = new JButton("Remove");
        remove.addActionListener(this);
        keywordValues.add(remove);
        JTextField textField = new JTextField();
        textField.setPreferredSize(entryDimension);
        textField.setMaximumSize(entryDimension);
        keywordValues.add(textField);
        JButton add = new JButton("Add");
        add.addActionListener(this);
        keywordValues.add(add);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Overidden equals method return true if object equals this, or if it of
     * the same class and has the same Tinker Keyword.
     */
    @Override
    public boolean equals(Object object) {
        if (this == object) {
            return true;
        } else if (object == null || getClass() != object.getClass()) {
            return false;
        }
        KeywordComponent other = (KeywordComponent) object;
        return keyword.equalsIgnoreCase(other.keyword);
    }

    /**
     * <p>
     * Getter for the field <code>keyword</code>.</p>
     *
     * @return a {@link java.lang.String} object.
     */
    public String getKeyword() {
        return keyword;
    }

    /**
     * <p>
     * getKeywordData</p>
     *
     * @param keywordData a {@link ffx.utilities.Keyword} object.
     */
    void getKeywordData(Keyword keywordData) {
        synchronized (this) {
            if (keywordData == null || !active) {
                return;
            }
            for (Component c : keywordValues) {
                if (c instanceof JTextField) {
                    JTextField tf = (JTextField) c;
                    if (!tf.getText().equals("")) {
                        keywordData.append(tf.getText());
                    }
                    break;
                } else if (c instanceof JComboBox) {
                    JComboBox cb = (JComboBox) c;
                    if (swingRepresentation == SwingRepresentation.COMBOBOX) {
                        String s = (String) cb.getSelectedItem();
                        if ("DEFAULT".equals(s)) {
                            logger.log(Level.WARNING, "Keyword should not be active: {0}", toString());
                            return;
                        }
                        keywordData.append(s);
                    } else {
                        int num = cb.getItemCount();
                        for (int i = 0; i < num; i++) {
                            String s = (String) cb.getItemAt(i);
                            keywordData.append(s);
                        }
                    }
                    break;
                } else if (c instanceof JCheckBox) {
                    JCheckBox cb = (JCheckBox) c;
                    String text = cb.getText();
                    if (cb.isSelected()) {
                        keywordData.append(text);
                    }
                }
            }
        }
    }

    /**
     * <p>
     * Getter for the field <code>keywordDescription</code>.</p>
     *
     * @return a {@link java.lang.String} object.
     */
    String getKeywordDescription() {
        return keywordDescription;
    }

    /**
     * <p>
     * Getter for the field <code>keywordGroup</code>.</p>
     *
     * @return a {@link java.lang.String} object.
     */
    String getKeywordGroup() {
        return keywordGroup;
    }

    /**
     * Returns a JPanel with a GridLayout LayoutManager that contains a Swing
     * representation of the Keyword and Modifiers in a single row.
     *
     * @return a {@link javax.swing.JPanel} object.
     */
    JPanel getKeywordGUI() {
        synchronized (this) {
            if (keywordGUI == null) {
                if (!init) {
                    initSwingComponents();
                }
                if (swingRepresentation == SwingRepresentation.MULTIPOLE) {
                    keywordGUI.add(keywordValues.get(0));
                } else {
                    keywordGUI = new JPanel(flowLayout);
                    for (Component c : keywordValues) {
                        keywordGUI.add(c);
                    }
                }
            }
            return keywordGUI;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int hashCode() {
        return Objects.hash(keyword.hashCode());
    }

    private void initSwingComponents() {
        if (swingRepresentation == SwingRepresentation.CHECKBOX) {
            checkBoxInit();
        } else if (swingRepresentation == SwingRepresentation.CHECKBOXES) {
            checkBoxesInit();
        } else if (swingRepresentation == SwingRepresentation.COMBOBOX) {
            comboBoxInit();
        } else if (swingRepresentation == SwingRepresentation.EDITCOMBOBOX) {
            editComboBoxInit();
        } else if (swingRepresentation == SwingRepresentation.TEXTFIELD) {
            textFieldInit();
        } else {
            return;
        }
        init = true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void insertUpdate(DocumentEvent evt) {
        isModified = true;
    }

    /**
     * <p>
     * isActive</p>
     *
     * @return a boolean.
     */
    public boolean isActive() {
        return active;
    }

    /**
     * Load a single line Keyword entry into this KeywordComponent. Keywords
     * that can be repeated multipule times are ComboBoxes are stored in
     * ComboBoxes.
     *
     * @param s A Keyword line, not including the Keyword itself.
     */
    @SuppressWarnings("unchecked")
    void loadKeywordEntry(String s) {
        synchronized (this) {
            if (!init) {
                initSwingComponents();
            }
            for (Component c : keywordValues) {
                if (c instanceof JCheckBox) {
                    JCheckBox cb = (JCheckBox) c;
                    if (s != null && s.equalsIgnoreCase(cb.getText())) {
                        cb.setSelected(true);
                    } else if (s == null) {
                        cb.setSelected(false);
                    }
                } else if (c instanceof JTextField) {
                    JTextField tf = (JTextField) c;
                    tf.setText(s);
                    break;
                } else if (c instanceof JComboBox && s != null) {
                    JComboBox<String> cb = (JComboBox<String>) c;
                    if (swingRepresentation == SwingRepresentation.EDITCOMBOBOX) {
                        cb.addItem(s);
                    } else {
                        cb.setSelectedItem(s);
                    }
                    break;
                }
            }
            active = toString() != null;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void mouseClicked(MouseEvent evt) {
        synchronized (this) {
            active = toString() != null;
            output.setText(keyword + ": " + keywordDescription);
            JViewport jsp = (JViewport) output.getParent();
            jsp.setViewPosition(new Point(20, 20));
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void mouseEntered(MouseEvent evt) {
        mouseClicked(evt);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void mouseExited(MouseEvent evt) {
        mouseClicked(evt);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void mousePressed(MouseEvent evt) {
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void mouseReleased(MouseEvent evt) {
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void removeUpdate(DocumentEvent evt) {
        isModified = true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void stateChanged(ChangeEvent evt) {
        isModified = true;
    }

    private void textFieldInit() {
        JLabel jl = new JLabel(keyword);
        jl.addMouseListener(this);
        jl.setPreferredSize(labelDimension);
        jl.setMaximumSize(labelDimension);
        keywordValues.add(jl);
        JTextField tf;
        tf = new JTextField("");
        tf.setColumns(25);
        tf.setActionCommand(keyword);
        tf.addMouseListener(this);
        tf.getDocument().addDocumentListener(this);
        if (keyword.equalsIgnoreCase("PARAMETERS")
                || keyword.equalsIgnoreCase("FORCEFIELD")) {
            tf.setEditable(false);
            EmptyBorder b = new EmptyBorder(1, 1, 1, 1);
            tf.setBorder(b);
            tf.setForeground(Color.BLUE);
        }
        tf.setPreferredSize(entryDimension);
        tf.setMaximumSize(entryDimension);
        keywordValues.add(tf);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Overridden toString methods facilitates Keyword output to a file.
     */
    @Override
    public String toString() {
        synchronized (this) {
            if (!active || !init) {
                return null;
            }
            StringBuilder s;
            // Torsion entries are long...
            // Some static layout variables
            String spaces = "                             ";
            if (!keyword.equalsIgnoreCase("TORSION")) {
                s = new StringBuilder(keyword
                        + spaces.substring(0, 18 - keyword.length()));
            } else {
                s = new StringBuilder(keyword);
            }
            for (Component c : keywordValues) {
                if (c instanceof JCheckBox) {
                    JCheckBox cb = (JCheckBox) c;
                    if (keywordValues.size() == 1) {
                        if (!cb.isSelected()) {
                            return null;
                        }
                    } else if (cb.getText().equalsIgnoreCase(keyword)) {
                        if (!cb.isSelected()) {
                            s = new StringBuilder();
                        }
                    } else {
                        if (cb.isSelected()) {
                            if (s.length() > 0) {
                                s.append("\n").append(keyword).append(" ").append(cb.getText());
                            } else {
                                s.append(keyword).append(" ").append(cb.getText());
                            }
                        }
                    }
                } else if (c instanceof JTextField) {
                    JTextField tf = (JTextField) c;
                    if (tf.getText().equals("")) {
                        return null;
                    }
                    String v = tf.getText();
                    if (v.length() < 8) {
                        s.append(v).append(spaces, 0, 8 - v.length());
                    } else {
                        s.append(v);
                    }
                    break;
                } else if (c instanceof JComboBox) {
                    JComboBox cb = (JComboBox) c;
                    if (swingRepresentation == SwingRepresentation.EDITCOMBOBOX) {
                        int count = cb.getItemCount();
                        if (count == 0) {
                            return null;
                        }
                        String[] entries = new String[count];
                        for (int i = 0; i < count; i++) {
                            entries[i] = (String) cb.getItemAt(i);
                        }
                        java.util.Arrays.sort(entries);
                        StringBuilder sb = new StringBuilder();
                        for (int i = 0; i < count; i++) {
                            sb.append(keyword).append(spaces.substring(0, 18 - keyword.length()));
                            sb.append(entries[i].toUpperCase());
                            if (i < count - 1) {
                                sb.append("\n");
                            }
                        }
                        s = sb;
                    } else {
                        String selection = (String) cb.getSelectedItem();
                        if (selection.startsWith("DEFAULT")) {
                            return null;
                        } else if (!selection.equalsIgnoreCase("PRESENT")) {
                            s.append(" ").append(selection);
                        }
                    }
                    break;
                }
            }
            if (s.length() == 0) {
                return null;
            }
            return s.toString();
        }
    }
}
