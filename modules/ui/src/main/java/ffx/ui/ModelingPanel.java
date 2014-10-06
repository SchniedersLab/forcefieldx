/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 */
package ffx.ui;

import java.awt.BorderLayout;
import java.awt.Button;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.prefs.Preferences;

import static java.lang.String.format;
import static java.util.Arrays.copyOfRange;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.JToggleButton;
import javax.swing.JToolBar;
import javax.swing.border.Border;
import javax.swing.border.EtchedBorder;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.apache.commons.io.FilenameUtils;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import ffx.potential.Utilities.FileType;
import ffx.potential.bonded.Residue;
import ffx.potential.parsers.SystemFilter;
import ffx.ui.commands.DTDResolver;
import ffx.utilities.Keyword;
import ffx.utilities.StringUtils;

/**
 * The ModelingPanel class encapsulates functionality needed to run FFX Modeling
 * Commands.
 *
 * @author Michael J. Schnieders
 *
 */
public class ModelingPanel extends JPanel implements ActionListener,
        MouseListener {

    private static final Logger logger = Logger.getLogger(ModelingPanel.class.getName());
    private static final Preferences prefs = Preferences.userNodeForPackage(ModelingPanel.class);
    private static final long serialVersionUID = 1L;

    private final MainPanel mainPanel;
    /**
     * Active System
     */
    private FFXSystem activeSystem = null;
    /**
     * File Type for the Active System
     */
    private FileType activeFileType = null;
    /**
     * Currently Selected Command
     */
    private String activeCommand = null;
    /**
     * File Types for this Command
     */
    private final ArrayList<FileType> commandFileTypes = new ArrayList<>();
    /**
     * Actions to take for this Command when it finishes
     */
    private String commandActions = "NONE";
    /**
     * Executing Commands
     */
    private final ArrayList<Thread> executingCommands = new ArrayList<>();
    /**
     * Log Settings
     */
    private final JComboBox<String> logSettings = new JComboBox<>();
    private String logString = null;
    /**
     * Commands for supported file types
     */
    private NodeList commandList;
    private JComboBox<String> xyzCommands;
    private JComboBox<String> intCommands;
    private JComboBox<String> arcCommands;
    private JComboBox<String> pdbCommands;
    private JComboBox<String> anyCommands;
    private JComboBox<String> currentCommandBox;
    /**
     * The CommandPanel holds the toolBar (north), splitPane (center) and
     * statusLabel (south).
     */
    private JPanel commandPanel;
    private JToolBar toolBar;
    /**
     * The splitPane holds the optionsTabbedPane (top) and descriptScrollPane
     * (bottom).
     */
    private JSplitPane splitPane;
    private JLabel statusLabel;
    private JTabbedPane optionsTabbedPane;
    private JScrollPane descriptScrollPane;
    private JTextArea descriptTextArea;
    private JCheckBoxMenuItem descriptCheckBox;
    private final ArrayList<JLabel> conditionals = new ArrayList<>();
    /**
     * Command input is formed in the commandTextArea, then exported to an input
     * file.
     */
    private JTextArea commandTextArea;
    /**
     * Nucleic Acid and Protein builder components.
     */
    private final JComboBox<String> acidComboBox = new JComboBox<>();
    private final JTextField acidTextField = new JTextField();
    private final JTextArea acidTextArea = new JTextArea();
    private JComboBox<String> conformationComboBox;
    private JScrollPane acidScrollPane = null;
    private JPanel aminoPanel = null;
    private JPanel nucleicPanel = null;
    /**
     * Reused FlowLayout.
     */
    private final FlowLayout flowLayout = new FlowLayout(FlowLayout.LEFT, 5, 5);
    /**
     * Reused BorderLayout.
     */
    private final BorderLayout borderLayout = new BorderLayout();
    /**
     * Reused EtchedBorder.
     */
    private final Border etchedBorder = BorderFactory.createEtchedBorder(EtchedBorder.RAISED);
    private final JTextField sizer = new JTextField(20);

    private JButton jbLaunch;
    private JButton jbStop;
    private FFXLauncher ffxLauncher = null;
    private Thread ffxThread = null;

    /**
     * Constructor
     *
     * @param mainPanel a {@link ffx.ui.MainPanel} object.
     */
    public ModelingPanel(MainPanel mainPanel) {
        super();
        this.mainPanel = mainPanel;
        initialize();
    }

    /**
     * {@inheritDoc}
     *
     * @param evt
     */
    @Override
    public void actionPerformed(ActionEvent evt) {
        synchronized (this) {
            String actionCommand = evt.getActionCommand();
            // A change to the selected TINKER Command
            switch (actionCommand) {
                case "FFXCommand":
                    JComboBox jcb = (JComboBox) toolBar.getComponentAtIndex(2);
                    String com = jcb.getSelectedItem().toString();
                    if (!com.equals(activeCommand)) {
                        activeCommand = com.toLowerCase();
                        loadCommand();
                    }
                    break;
                case "LogSettings":
                    // A change to the Log Settings.
                    loadLogSettings();
                    statusLabel.setText("  " + createCommandInput());
                    break;
                case "Launch":
                    // Launch the selected Force Field X command.
                    runScript();
                    break;
                case "NUCLEIC":
                case "PROTEIN":
                    // Editor functions for the Protein and Nucleic Acid Builders
                    builderCommandEvent(evt);
                    break;
                case "Conditional":
                    // Some command options are conditional on other input.
                    conditionalCommandEvent(evt);
                    break;
                case "End":
                    // End the currently executing command.
                    setEnd();
                    break;
                case "Delete":
                    // Delete log files.
                    deleteLogs();
                    break;
                case "Description":
                    // Allow command descriptions to be hidden.
                    JCheckBoxMenuItem box = (JCheckBoxMenuItem) evt.getSource();
                    setDivider(box.isSelected());
                    break;
                default:
                    logger.log(Level.WARNING, "ModelingPanel ActionCommand not recognized: {0}", evt);
                    break;
            }
        }
    }

    /*
     * This handles Protein and Nucleic Acid builder events.
     */
    private void builderCommandEvent(ActionEvent evt) {
        JButton button = (JButton) evt.getSource();
        String arg = evt.getActionCommand();
        int index = acidComboBox.getSelectedIndex();
        String selected = (String) acidComboBox.getItemAt(index);
        if ("Remove".equals(button.getText())) {
            // Remove one entry
            if (acidComboBox.getItemCount() > 0) {
                acidComboBox.removeItemAt(index);
                index--;
            }
        } else if ("Edit".equals(button.getText())) {
            String entry = new String(acidTextField.getText());
            // Allow editing - should add more input validation here
            if (!entry.equals("")) {
                String s[] = entry.trim().split(" +");
                String newResidue = s[0].toUpperCase();
                if ("NUCLEIC".equals(arg)) {
                    // Residue.NA3Set.contains(newResidue);
                    try {
                        Residue.NA3.valueOf(newResidue);
                        acidComboBox.removeItemAt(index);
                        acidComboBox.insertItemAt("" + index + " " + entry,
                                index);
                    } catch (Exception e) {
                    }
                } else {
                    try {
                        Residue.AA3.valueOf(newResidue);
                        acidComboBox.removeItemAt(index);
                        acidComboBox.insertItemAt("" + index + " " + entry,
                                index);
                    } catch (Exception e) {
                    }
                }
            }
        } else if ("Reset".equals(button.getText())) {
            // Remove all entries
            acidComboBox.removeAllItems();
            acidTextArea.setText("");
        } else {
            // A base/residue button was selected
            String newResidue = button.getText();
            if ("PROTEIN".equals(arg)) {
                String c = (String) conformationComboBox.getSelectedItem();
                if (!c.toUpperCase().startsWith("DEFAULT")) {
                    c = c.substring(c.indexOf("[") + 1, c.indexOf("]"));
                    newResidue = newResidue + " " + c;
                }
                acidComboBox.insertItemAt("" + index + " " + newResidue,
                        index + 1);
                index++;
            } else {
                if (!newResidue.equalsIgnoreCase("MOL")) {
                    acidComboBox.insertItemAt("" + index + " " + newResidue,
                            index + 1);
                    index++;
                } else if (!selected.equalsIgnoreCase("MOL")) {
                    acidComboBox.insertItemAt("" + index + " " + newResidue,
                            index + 1);
                    index++;
                }
            }
        }
        // Create the condensed sequence view.
        StringBuilder sequence = new StringBuilder();
        for (int i = 0; i < acidComboBox.getItemCount(); i++) {
            String s[] = ((String) acidComboBox.getItemAt(i)).trim().toUpperCase().split(" +");
            if (s.length > 1) {
                if (s[1].equalsIgnoreCase("MOL")) {
                    sequence.append(s[1]).append("\n");
                } else {
                    sequence.append(s[1]).append(" ");
                }
            }
        }
        // Renumber the sequence.
        acidTextArea.setText(sequence.toString());
        for (int i = 0; i < acidComboBox.getItemCount(); i++) {
            String s = (String) acidComboBox.getItemAt(i);
            s = s.substring(s.indexOf(" "), s.length()).trim();
            acidComboBox.removeItemAt(i);
            acidComboBox.insertItemAt("" + (i + 1) + " " + s, i);
        }
        // Set the selected entry and fill in the edit textField.
        if (index < 0) {
            index = 0;
        }
        if (index > acidComboBox.getItemCount() - 1) {
            index = acidComboBox.getItemCount() - 1;
        }
        acidComboBox.setSelectedIndex(index);
        String s = (String) acidComboBox.getItemAt(index);
        if (s != null) {
            acidTextField.setText(s.substring(s.indexOf(" "), s.length()).trim());
        } else {
            acidTextField.setText("");
        }
    }

    /**
     * This handles conditional command option input.
     *
     * @param evt ActionEvent
     */
    private void conditionalCommandEvent(ActionEvent evt) {
        Object source = evt.getSource();
        if (source instanceof JRadioButton) {
            JRadioButton jrb = (JRadioButton) source;
            String selection = jrb.getText().toLowerCase();
            for (JLabel label : conditionals) {
                JTextField jtf = (JTextField) label.getLabelFor();
                String cupon = label.getName().toLowerCase();
                if (cupon.contains(selection) && jrb.isSelected()) {
                    label.setEnabled(true);
                    jtf.setEnabled(true);
                } else {
                    label.setEnabled(false);
                    jtf.setEnabled(false);
                }
            }
        } else if (source instanceof JCheckBox) {
            JCheckBox jcb = (JCheckBox) source;
            String selection = jcb.getText().toLowerCase();
            for (JLabel label : conditionals) {
                String cupon = label.getName().toLowerCase();
                JTextField jtf = (JTextField) label.getLabelFor();
                if (cupon.contains(selection) && jcb.isSelected()) {
                    label.setEnabled(true);
                    jtf.setEnabled(true);
                } else {
                    label.setEnabled(false);
                    jtf.setEnabled(false);
                }
            }
        }
        statusLabel.setText("  " + createCommandInput());
    }

    /**
     * Create a string representing the modeling command to execute.
     *
     * @return
     */
    private String createCommandInput() {
        StringBuilder commandLineParams = new StringBuilder(activeCommand + " ");
        // Now append command line input to a TextArea, one option per line.
        // This TextArea gets dumped to an input file.
        commandTextArea.setText("");
        int numparams = optionsTabbedPane.getTabCount();
        for (int i = 0; i < numparams; i++) {
            // A few cases require that a newLine not be generated between
            // options.
            boolean newLine = true;
            // The optionString will collect the parameters for this Option,
            // then append them to the CommandTextArea.
            StringBuilder optionString = new StringBuilder();
            JPanel optionPanel = (JPanel) optionsTabbedPane.getComponentAt(i);
            int numOptions = optionPanel.getComponentCount();
            String title = optionsTabbedPane.getTitleAt(i);
            if (title.equalsIgnoreCase("Sequence")) {
                for (int k = 0; k < acidComboBox.getItemCount(); k++) {
                    if (k != 0) {
                        optionString.append("\n");
                    }
                    String s = (String) acidComboBox.getItemAt(k);
                    s = s.substring(s.indexOf(" "), s.length()).trim();
                    optionString.append(s);
                }
                // Need an extra newline for Nucleic
                if (activeCommand.equalsIgnoreCase("NUCLEIC")) {
                    optionString.append("\n");
                }
            } else {
                JPanel valuePanel = (JPanel) optionPanel.getComponent(numOptions - 1);
                int numValues = valuePanel.getComponentCount();
                for (int j = 0; j < numValues; j++) {
                    Component value = valuePanel.getComponent(j);
                    if (value instanceof JCheckBox) {
                        JCheckBox jcbox = (JCheckBox) value;
                        if (jcbox.isSelected()) {
                            optionString.append("-");
                            optionString.append(jcbox.getName());
                            optionString.append(" ");
                            optionString.append(jcbox.getText());
                        }
                    } else if (value instanceof JTextField) {
                        JTextField jtfield = (JTextField) value;
                        optionString.append("-");
                        optionString.append(jtfield.getName());
                        optionString.append(" ");
                        optionString.append(jtfield.getText());
                    } else if (value instanceof JComboBox) {
                        JComboBox jcb = (JComboBox) value;
                        Object object = jcb.getSelectedItem();
                        if (object instanceof FFXSystem) {
                            FFXSystem system = (FFXSystem) object;
                            File file = system.getFile();
                            if (file != null) {
                                String absolutePath = file.getAbsolutePath();
                                if (absolutePath.endsWith("xyz")) {
                                    absolutePath = absolutePath + "_1";
                                }
                                optionString.append(absolutePath);
                            }
                        }
                    } else if (value instanceof JRadioButton) {
                        JRadioButton jrbutton = (JRadioButton) value;
                        if (jrbutton.isSelected()) {
                            if (!jrbutton.getText().equalsIgnoreCase("NONE")) {
                                optionString.append("-");
                                optionString.append(jrbutton.getName());
                                optionString.append(" ");
                                optionString.append(jrbutton.getText());
                            }
                            if (title.equalsIgnoreCase("C-CAP")) {
                                optionString.append("\n");
                            }
                        }
                    }
                }
                // Handle Conditional Options
                if (optionPanel.getComponentCount() == 3) {
                    valuePanel = (JPanel) optionPanel.getComponent(1);
                    // JLabel conditionalLabel = (JLabel)
                    // valuePanel.getComponent(0);
                    JTextField jtf = (JTextField) valuePanel.getComponent(1);
                    if (jtf.isEnabled()) {
                        String conditionalInput = jtf.getText();
                        // Post-Process the Input into Atom Pairs
                        String postProcess = jtf.getName();
                        if (postProcess != null
                                && postProcess.equalsIgnoreCase("ATOMPAIRS")) {
                            String tokens[] = conditionalInput.split(" +");
                            StringBuilder atomPairs = new StringBuilder();
                            int atomNumber = 0;
                            for (String token : tokens) {
                                atomPairs.append(token);
                                if (atomNumber++ % 2 == 0) {
                                    atomPairs.append(" ");
                                } else {
                                    atomPairs.append("\n");
                                }
                            }
                            conditionalInput = atomPairs.toString();
                        }
                        // Append a newline to "enter" the option string.
                        // Append "conditional" input.
                        optionString.append("\n").append(conditionalInput);
                    }
                }
            }
            if (optionString.length() > 0) {
                commandTextArea.append(optionString.toString());
                if (newLine) {
                    commandTextArea.append("\n");
                }
            }
        }
        String commandInput = commandTextArea.getText();
        if (commandInput != null && !commandInput.trim().equalsIgnoreCase("")) {
            commandLineParams.append(commandInput);
        }

        // The final token on the command line is the structure file name, except
        // for protein and nucleic.
        if (!activeCommand.equalsIgnoreCase("Protein")
                && !activeCommand.equalsIgnoreCase("Nucleic")) {
            File file = activeSystem.getFile();
            if (file != null) {
                String name = file.getName();
                commandLineParams.append(name);
                commandLineParams.append(" ");
            } else {
                return null;
            }
        }

        return commandLineParams.toString();
    }

    /**
     * Launch the TINKER command specified by the ModelingPanel
     *
     * @return a {@link ffx.ui.FFXExec} object.
     */
    public FFXExec executeCommand() {
        FFXSystem s = mainPanel.getHierarchy().getActive();
        String dir = MainPanel.getPWD().getAbsolutePath();
        if (s != null) {
            File f = s.getFile();
            if (f != null) {
                dir = f.getParent();
            }
        }
        return launch(statusLabel.getText(), dir);
    }

    private JPanel getAminoAcidPanel() {
        if (aminoPanel != null) {
            return aminoPanel;
        }
        JPanel buttonPanel = new JPanel(new GridLayout(4, 5, 2, 2));
        Residue.AA3 a[] = Residue.AA3.values();
        for (int i = 0; i < 20; i++) {
            JButton button = new JButton(a[i].name());
            button.setActionCommand("PROTEIN");
            button.addActionListener(this);
            buttonPanel.add(button);
        }
        buttonPanel.setMaximumSize(buttonPanel.getPreferredSize());
        aminoPanel = new JPanel();
        aminoPanel.setLayout(new BoxLayout(aminoPanel, BoxLayout.Y_AXIS));
        aminoPanel.add(buttonPanel);
        conformationComboBox = new JComboBox<>(Residue.Ramachandran);
        conformationComboBox.setFont(Font.decode("Monospaced"));
        conformationComboBox.setMaximumSize(buttonPanel.getPreferredSize());
        aminoPanel.add(conformationComboBox);
        return aminoPanel;
    }

    /**
     * <p>
     * getAvailableCommands</p>
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<String> getAvailableCommands() {
        ArrayList<String> availableCommands = new ArrayList<>();
        for (int i = 0; i < currentCommandBox.getItemCount(); i++) {
            availableCommands.add((String) currentCommandBox.getItemAt(i));
        }
        return availableCommands;
    }

    /**
     **********************************************************************
     *
     * @return a {@link java.lang.String} object.
     */
    // Modeling Command Configuration
    public String getCommand() {
        return activeCommand;
    }

    private String getLogString(File currentLog) {
        String currentMode = (String) logSettings.getSelectedItem();
        if (currentMode.startsWith("Create")) {
            currentLog = SystemFilter.version(currentLog);
            return " > \"" + currentLog.getAbsolutePath() + "\"";
        } else if (currentMode.startsWith("Append")) {
            return " >> \"" + currentLog.getAbsolutePath() + "\"";
        } else {
            return " > \"" + currentLog.getAbsolutePath() + "\"";
        }
    }

    /**
     * Get an ArrayList of executing jobs
     *
     * @return an ArrayList containing Thread objects
     */
    public ArrayList<Thread> getModelingJobs() {
        return executingCommands;
    }

    private JPanel getNucleicAcidPanel() {
        if (nucleicPanel != null) {
            return nucleicPanel;
        }
        nucleicPanel = new JPanel(new GridLayout(3, 4, 2, 2));
        Residue.NA3 a[] = Residue.NA3.values();
        for (int i = 0; i < 8; i++) {
            JButton button = new JButton(a[i].name());
            button.setActionCommand("NUCLEIC");
            button.addActionListener(this);
            nucleicPanel.add(button);
        }
        JButton button = new JButton("MOL");
        button.setActionCommand("NUCLEIC");
        button.addActionListener(this);
        nucleicPanel.add(button);
        nucleicPanel.add(Box.createHorizontalBox());
        nucleicPanel.add(Box.createHorizontalBox());
        nucleicPanel.add(Box.createHorizontalBox());
        nucleicPanel.setMaximumSize(nucleicPanel.getPreferredSize());
        return nucleicPanel;
    }

    private void initCommandComboBox(JComboBox commands) {
        commands.setActionCommand("FFXCommand");
        commands.setMaximumSize(xyzCommands.getPreferredSize());
        commands.setEditable(false);
        commands.setToolTipText("Select a Modeling Command");
        commands.setSelectedIndex(0);
        commands.addActionListener(this);
    }

    private void initialize() {
        // Command Description
        descriptTextArea = new JTextArea();
        descriptTextArea.setEditable(false);
        descriptTextArea.setLineWrap(true);
        descriptTextArea.setWrapStyleWord(true);
        descriptTextArea.setDoubleBuffered(true);
        Insets insets = descriptTextArea.getInsets();
        insets.set(5, 5, 5, 5);
        descriptTextArea.setMargin(insets);
        descriptScrollPane = new JScrollPane(descriptTextArea,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        descriptScrollPane.setBorder(etchedBorder);
        // Command Input
        commandTextArea = new JTextArea();
        commandTextArea.setEditable(false);
        commandTextArea.setLineWrap(true);
        commandTextArea.setWrapStyleWord(true);
        commandTextArea.setDoubleBuffered(true);
        commandTextArea.setMargin(insets);
        // Command Options
        optionsTabbedPane = new JTabbedPane();
        statusLabel = new JLabel();
        statusLabel.setBorder(etchedBorder);
        statusLabel.setToolTipText("  Modeling command that will be executed");
        commandPanel = new JPanel(flowLayout);
        commandPanel.add(optionsTabbedPane);
        splitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, commandPanel,
                descriptScrollPane);
        splitPane.setContinuousLayout(true);
        splitPane.setResizeWeight(1.0d);
        splitPane.setOneTouchExpandable(true);
        setLayout(new BorderLayout());
        add(splitPane, BorderLayout.CENTER);
        add(statusLabel, BorderLayout.SOUTH);
        // Initialize the Amino/Nucleic Acid ComboBox.
        acidComboBox.setEditable(false);
        acidComboBox.setMaximumSize(sizer.getPreferredSize());
        acidComboBox.setPreferredSize(sizer.getPreferredSize());
        acidComboBox.setMinimumSize(sizer.getPreferredSize());
        acidComboBox.setFont(Font.decode("Monospaced"));
        acidTextField.setMaximumSize(sizer.getPreferredSize());
        acidTextField.setMinimumSize(sizer.getPreferredSize());
        acidTextField.setPreferredSize(sizer.getPreferredSize());
        acidTextArea.setEditable(false);
        acidTextArea.setWrapStyleWord(true);
        acidTextArea.setLineWrap(true);
        acidTextArea.setFont(Font.decode("Monospaced"));
        acidScrollPane = new JScrollPane(acidTextArea);
        Dimension d = new Dimension(300, 400);
        acidScrollPane.setPreferredSize(d);
        acidScrollPane.setMaximumSize(d);
        acidScrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
        // Load the FFX commands.xml file that defines FFX commands.
        try {
            DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
            DocumentBuilder db = dbf.newDocumentBuilder();
            db.setEntityResolver(new DTDResolver());
            URL comURL = getClass().getClassLoader().getResource(
                    "ffx/ui/commands/commands.xml");
            Document doc = db.parse(comURL.openStream());
            NodeList nodelist = doc.getChildNodes();
            Node commandroot = null;
            for (int i = 0; i < nodelist.getLength(); i++) {
                commandroot = nodelist.item(i);
                if (commandroot.getNodeName().equals("FFXCommands")
                        && commandroot instanceof Element) {
                    break;
                }
            }
            if (commandroot == null || !(commandroot instanceof Element)) {
                commandList = null;
            }
            commandList = ((Element) commandroot).getElementsByTagName("Command");
        } catch (ParserConfigurationException | SAXException | IOException e) {
            System.err.println(e);
        } finally {
            if (commandList == null) {
                System.out.println("Force Field X commands.xml could not be parsed.");
                logger.severe("Force Field X will exit.");
                System.exit(-1);
            }
        }
        // Create a ComboBox with commands specific to each type of coordinate
        // file.
        xyzCommands = new JComboBox<>();
        intCommands = new JComboBox<>();
        arcCommands = new JComboBox<>();
        pdbCommands = new JComboBox<>();
        anyCommands = new JComboBox<>();
        Element command;
        String name;
        int numcommands = commandList.getLength();
        for (int i = 0; i < numcommands; i++) {
            command = (Element) commandList.item(i);
            name = command.getAttribute("name");
            String temp = command.getAttribute("fileType");
            if (temp.contains("ANY")) {
                temp = "XYZ INT ARC PDB";
                anyCommands.addItem(name);
            }
            String[] types = temp.split(" +");
            for (String type : types) {
                if (type.contains("XYZ")) {
                    xyzCommands.addItem(name);
                }
                if (type.contains("INT")) {
                    intCommands.addItem(name);
                }
                if (type.contains("ARC")) {
                    arcCommands.addItem(name);
                }
                if (type.contains("PDB")) {
                    pdbCommands.addItem(name);
                }
            }
        }
        initCommandComboBox(xyzCommands);
        initCommandComboBox(intCommands);
        initCommandComboBox(arcCommands);
        initCommandComboBox(pdbCommands);
        initCommandComboBox(anyCommands);
        currentCommandBox = anyCommands;
        activeCommand = (String) anyCommands.getSelectedItem();
        // Load the default Command.
        loadCommand();
        // Load the default Log File Settings.
        logSettings.setActionCommand("LogSettings");
        loadLogSettings();

        // Create the Toolbar.
        toolBar = new JToolBar("Modeling Commands", JToolBar.HORIZONTAL);
        toolBar.setLayout(new FlowLayout(FlowLayout.LEFT));
        jbLaunch = new JButton(new ImageIcon(getClass().getClassLoader().getResource("ffx/ui/icons/cog_go.png")));
        jbLaunch.setActionCommand("Launch");
        jbLaunch.setToolTipText("Launch the Force Field X Command");
        jbLaunch.addActionListener(this);
        insets.set(2, 2, 2, 2);
        jbLaunch.setMargin(insets);
        toolBar.add(jbLaunch);
        jbStop = new JButton(new ImageIcon(getClass().getClassLoader().getResource("ffx/ui/icons/stop.png")));
        jbStop.setActionCommand("End");
        jbStop.setToolTipText("Terminate the Current Force Field X Command");
        jbStop.addActionListener(this);
        jbStop.setMargin(insets);
        jbStop.setEnabled(false);
        toolBar.add(jbStop);
        toolBar.addSeparator();
        toolBar.add(anyCommands);
        currentCommandBox = anyCommands;
        toolBar.addSeparator();
        /*
         toolBar.add(logSettings);
         JButton jbdelete = new JButton(new ImageIcon(getClass().getClassLoader().getResource("ffx/ui/icons/page_delete.png")));
         jbdelete.setActionCommand("Delete");
         jbdelete.setToolTipText("Delete Log Files");
         jbdelete.addActionListener(this);
         jbdelete.setMargin(insets);
         toolBar.add(jbdelete);
         toolBar.addSeparator();
         */
        ImageIcon icinfo = new ImageIcon(getClass().getClassLoader().getResource("ffx/ui/icons/information.png"));
        descriptCheckBox = new JCheckBoxMenuItem(icinfo);
        descriptCheckBox.addActionListener(this);
        descriptCheckBox.setActionCommand("Description");
        descriptCheckBox.setToolTipText("Show/Hide Modeling Command Descriptions");
        descriptCheckBox.setMargin(insets);
        toolBar.add(descriptCheckBox);
        toolBar.add(new JLabel(""));
        toolBar.setBorderPainted(false);
        toolBar.setFloatable(false);
        toolBar.setRollover(true);
        add(toolBar, BorderLayout.NORTH);
        // Load ModelingPanel preferences.
        Preferences prefs = Preferences.userNodeForPackage(ffx.ui.ModelingPanel.class);
        descriptCheckBox.setSelected(!prefs.getBoolean("JobPanel_description",
                true));
        descriptCheckBox.doClick();
    }

    private void runScript() {
        String name = activeCommand;
        name = name.replace('.', File.separatorChar);
        ClassLoader loader = getClass().getClassLoader();
        URL embeddedScript = loader.getResource("ffx/scripts/" + name + ".ffx");
        if (embeddedScript == null) {
            embeddedScript = loader.getResource("ffx/scripts/" + name + ".groovy");
        }
        File scriptFile = null;
        if (embeddedScript != null) {
            try {
                scriptFile = new File(
                    StringUtils.copyInputStreamToTmpFile(
                                embeddedScript.openStream(), name, ".ffx"));
            } catch (IOException e) {
                logger.log(Level.WARNING, "Exception extracting embedded script {0}\n{1}",
                        new Object[]{embeddedScript.toString(), e.toString()});
            }
        }
        if (scriptFile != null && scriptFile.exists()) {
            String argLine = statusLabel.getText().replace('\n', ' ');
            String args[] = argLine.trim().split(" +");
            // Remove the command (first token) and system name (last token).
            args = copyOfRange(args, 1, args.length - 1);
            List<String> argList = Arrays.asList(args);
            ffxLauncher = new FFXLauncher(argList, scriptFile);
            ffxThread = new Thread(ffxLauncher);
            ffxThread.setName(statusLabel.getText());
            ffxThread.setPriority(Thread.MAX_PRIORITY);
            ffxThread.start();
        } else {
            logger.warning(format("%s was not found.", name));
        }
    }

    public void enableLaunch(boolean enable) {
        jbLaunch.setEnabled(enable);
        jbStop.setEnabled(!enable);
    }

    private class FFXLauncher implements Runnable {

        List<String> argList;
        File scriptFile;

        public FFXLauncher(List<String> argList, File scriptFile) {
            this.argList = argList;
            this.scriptFile = scriptFile;
        }

        @Override
        public void run() {
            mainPanel.getModelingShell().setArgList(argList);
            mainPanel.open(scriptFile, null);
        }
    }

    /**
     * Launch the active command on the active system in the specified
     * directory.
     *
     * @param command The command to be excuted.
     * @param dir The directory to execute the command in.
     * @return a {@link ffx.ui.FFXExec} object.
     */
    public FFXExec launch(String command, String dir) {
        logger.log(Level.INFO, "Command: {0}\nDirectory: {1}", new Object[]{command, dir});
        synchronized (this) {
            // Check that the TINKER *.exe exists in TINKER/bin
            String path = MainPanel.ffxDir.getAbsolutePath();
            File exe = new File(path + File.separator
                    + activeCommand.toLowerCase());
            if (!exe.exists()) {
                exe = new File(exe.getAbsolutePath() + ".exe");
                if (!exe.exists()) {
                    String message = "The "
                            + activeCommand
                            + " executable was not found in "
                            + path
                            + ". Please use the 'Set TINKER...' dialog to change the TINKER directory.";
                    JOptionPane.showMessageDialog(null, message,
                            "Could not launch " + activeCommand,
                            JOptionPane.ERROR_MESSAGE);
                    return null;
                }
            }
            // Check that the directory to execute the command in is valid
            File dirf = new File(dir);
            if (!dirf.exists()) {
                logger.log(Level.WARNING, "Directory doesn''t exist: {0}", dirf.getAbsolutePath());
                return null;
            }
            // Check if we need a key file
            if (!commandFileTypes.contains(FileType.ANY)) {
                if (activeSystem == null) {
                    return null;
                }
                activeFileType = FileType.XYZ;
                // Check that the TINKER command executes on this file type
                if (!commandFileTypes.contains(activeFileType)) {
                    String message = activeCommand.toUpperCase()
                            + " does not execute on " + activeFileType
                            + " files.";
                    JOptionPane.showMessageDialog(null, message,
                            "Could not launch " + activeCommand,
                            JOptionPane.ERROR_MESSAGE);
                    return null;
                }
                // Check that a key file exists or prompt to create one.
                if (activeSystem.getKeyFile() == null) {
                    mainPanel.createKeyFile(activeSystem);
                    // Give up if the key file is null.
                    if (activeSystem.getKeyFile() == null) {
                        return null;
                    }
                }
            } else {
                // Determine names to use for the output of Protein/Nucleic
                command = createCommandInput();
                String structureName = commandTextArea.getText().trim();
                if (!structureName.equalsIgnoreCase("")) {
                    structureName = (structureName.split("\n"))[0];
                    if (structureName != null) {
                        structureName = structureName.trim();
                        int dot = structureName.lastIndexOf(".");
                        if (dot > 0) {
                            structureName = structureName.substring(0, dot);
                        }
                    }
                }
                // If the above fails, just use the name of the executable
                // (protein or nulceic)
                if (structureName == null) {
                    structureName = activeCommand.toLowerCase();
                }
                File file = new File(dir + File.separator + structureName
                        + ".xyz");
                file = SystemFilter.version(file);
                activeSystem = new FFXSystem(file, null, Keyword.loadProperties(file));
                File logFile = new File(file.getParent() + File.separator
                        + structureName + ".log");
                activeSystem.setLogFile(logFile);
                loadLogSettings();
                activeFileType = FileType.ANY;
                // Need to have a parameter file chosen.
                mainPanel.openKey(activeSystem, true);
                if (activeSystem.getKeyFile() == null) {
                    return null;
                }
            }
            // Decide on a Log file
            if (((String) logSettings.getSelectedItem()).startsWith("Create")) {
                File newLog = SystemFilter.version(activeSystem.getLogFile());
                activeSystem.setLogFile(newLog);
            }
            String logName = activeSystem.getLogFile().getAbsolutePath();
            // Determine the command string
            command = createCommandInput();
            // If a new structure file will be created, determine what the name
            // will be.
            File newFile = null;
            if (commandActions.toUpperCase().contains("LOAD")) {
                File oldFile = activeSystem.getFile();
                if (commandActions.toUpperCase().contains("LOADXYZ")) {
                    String fileName = oldFile.getAbsolutePath();
                    int dot = fileName.lastIndexOf(".");
                    if (dot > 0) {
                        fileName = fileName.substring(0, dot) + ".xyz";
                    }
                    oldFile = new File(fileName);
                } else if (commandActions.toUpperCase().contains("LOADINT")) {
                    String fileName = oldFile.getAbsolutePath();
                    int dot = fileName.lastIndexOf(".");
                    if (dot > 0) {
                        fileName = fileName.substring(0, dot) + ".int";
                    }
                    oldFile = new File(fileName);
                } else if (commandActions.toUpperCase().contains("LOADPDB")) {
                    String fileName = oldFile.getAbsolutePath();
                    int dot = fileName.lastIndexOf(".");
                    if (dot > 0) {
                        fileName = fileName.substring(0, dot) + ".pdb";
                    }
                    oldFile = new File(fileName);
                }
                newFile = SystemFilter.version(oldFile);
            }
            // Save any changes that have been made to the key file
            mainPanel.getKeywordPanel().saveChanges();
            // Remove any TINKER *.END files
            removeEnd();
            // Create the input file
            String commandInput = commandTextArea.getText();
            if (commandInput != null
                    && !commandInput.trim().equalsIgnoreCase("")) {
                File inputFile = new File(dir + File.separator
                        + activeCommand.toLowerCase() + ".in");
                inputFile.deleteOnExit();
                try {
                    FileWriter fw = new FileWriter(inputFile);
                    fw.write(commandInput);
                    fw.close();
                } catch (Exception e) {
                    logger.info(e.toString());
                    return null;
                }
            }
            // If the job progressively modifies coordinates, open a copy of it.
            boolean openOnto = false;
            if (commandActions.toUpperCase().contains("CONNECT")) {
                // If a version file is created, open it onto the structure used
                // to
                // display the job.
                if (newFile != null) {
                    openOnto = true;
                }
                mainPanel.open(activeSystem.getFile(), activeCommand);
                try {
                    while (mainPanel.isOpening()) {
                        wait(10);
                    }
                } catch (Exception e) {
                    logger.info(e.toString());
                    return null;
                }
                activeSystem = mainPanel.getHierarchy().getActive();
            }
            // Finally, create and execute the command in a new thread.
            FFXExec tinkerExec = new FFXExec(activeSystem, logName, command,
                    dir, mainPanel, newFile, openOnto);
            Thread tinkerThread = new Thread(tinkerExec);
            tinkerThread.setPriority(Thread.NORM_PRIORITY);
            tinkerThread.setName(logName);
            // If the job progressively modifies coordinates, connect to it.
            if (commandActions.toUpperCase().contains("CONNECT")) {
                mainPanel.connectToTINKER(activeSystem, tinkerThread);
            } else {
                tinkerThread.start();
            }
            // If some action should be taken when the job finishes,
            // add it to the Modeling Jobs Vector
            if (!commandActions.equalsIgnoreCase("NONE")) {
                executingCommands.add(tinkerThread);
                // mainPanel.getLogPanel().refreshStatus();
            }
            return tinkerExec;
        }
    }

    /**
     * Load the active system into the JobPanel. This should be called whenever
     * the active system changes.
     *
     * @param active a {@link ffx.ui.FFXSystem} object.
     */
    public void loadActive(FFXSystem active) {
        synchronized (this) {
            activeSystem = active;
            FileType fileType;
            // No Open Molecules
            if (activeSystem == null || activeSystem.isClosing()) {
                currentCommandBox = anyCommands;
                activeCommand = (String) anyCommands.getSelectedItem();
                statusLabel.setText("  ");
                fileType = FileType.ANY;
            } else {
                fileType = FileType.XYZ;
            }
            if (fileType != activeFileType) {
                activeFileType = fileType;
                toolBar.remove(2);
                if (activeFileType == FileType.XYZ) {
                    toolBar.add(xyzCommands);
                    currentCommandBox = xyzCommands;
                } else if (activeFileType == FileType.INT) {
                    toolBar.add(intCommands);
                    currentCommandBox = intCommands;
                } else if (activeFileType == FileType.ARC) {
                    toolBar.add(arcCommands);
                    currentCommandBox = arcCommands;
                } else if (activeFileType == FileType.PDB) {
                    toolBar.add(pdbCommands);
                    currentCommandBox = pdbCommands;
                } else {
                    toolBar.add(anyCommands);
                    currentCommandBox = anyCommands;
                }
                toolBar.add(currentCommandBox, 2);
                toolBar.validate();
                toolBar.repaint();
                activeCommand = (String) currentCommandBox.getSelectedItem();
            }
            loadCommand();
        }
    }

    private void loadCommand() {
        synchronized (this) {
            // Force Field X Command
            Element command;
            // Command Options
            NodeList options;
            Element option;
            // Option Values
            NodeList values;
            Element value;
            // Options may be Conditional based on previous Option values (not
            // always supplied)
            NodeList conditionalList;
            Element conditional;
            // JobPanel GUI Components that change based on command
            JPanel optionPanel;
            // Clear the previous components
            commandPanel.removeAll();
            optionsTabbedPane.removeAll();
            conditionals.clear();
            String currentCommand = (String) currentCommandBox.getSelectedItem();
            if (currentCommand == null) {
                commandPanel.validate();
                commandPanel.repaint();
                return;
            }
            command = null;
            for (int i = 0; i < commandList.getLength(); i++) {
                command = (Element) commandList.item(i);
                String name = command.getAttribute("name");
                if (name.equalsIgnoreCase(currentCommand)) {
                    break;
                }
            }
            int div = splitPane.getDividerLocation();
            descriptTextArea.setText(currentCommand.toUpperCase() + ": "
                    + command.getAttribute("description"));
            splitPane.setBottomComponent(descriptScrollPane);
            splitPane.setDividerLocation(div);
            // The "action" tells Force Field X what to do when the
            // command finishes
            commandActions = command.getAttribute("action").trim();
            // The "fileType" specifies what file types this command can execute
            // on
            String string = command.getAttribute("fileType").trim();
            String[] types = string.split(" +");
            commandFileTypes.clear();
            for (String type : types) {
                if (type.contains("XYZ")) {
                    commandFileTypes.add(FileType.XYZ);
                }
                if (type.contains("INT")) {
                    commandFileTypes.add(FileType.INT);
                }
                if (type.contains("ARC")) {
                    commandFileTypes.add(FileType.ARC);
                }
                if (type.contains("PDB")) {
                    commandFileTypes.add(FileType.PDB);
                }
                if (type.contains("ANY")) {
                    commandFileTypes.add(FileType.ANY);
                }
            }
            // Determine what options are available for this command
            options = command.getElementsByTagName("Option");
            int length = options.getLength();
            for (int i = 0; i < length; i++) {
                // This Option will be enabled (isEnabled = true) unless a
                // Conditional disables it
                boolean isEnabled = true;
                option = (Element) options.item(i);
                conditionalList = option.getElementsByTagName("Conditional");
                /*
                 * Currently, there can only be 0 or 1 Conditionals per Option
                 * There are three types of Conditionals implemented. 1.)
                 * Conditional on a previous Option, this option may be
                 * available 2.) Conditional on input for this option, a
                 * sub-option may be available 3.) Conditional on the presence
                 * of keywords, this option may be available
                 */
                if (conditionalList != null) {
                    conditional = (Element) conditionalList.item(0);
                } else {
                    conditional = null;
                }
                // Get the command line flag
                String flag = option.getAttribute("flag").trim();
                // Get the description
                String optionDescript = option.getAttribute("description");
                JTextArea optionTextArea = new JTextArea("  "
                        + optionDescript.trim());
                optionTextArea.setEditable(false);
                optionTextArea.setLineWrap(true);
                optionTextArea.setWrapStyleWord(true);
                optionTextArea.setBorder(etchedBorder);
                // Get the default for this Option (if one exists)
                String defaultOption = option.getAttribute("default");
                // Option Panel
                optionPanel = new JPanel(new BorderLayout());
                optionPanel.add(optionTextArea, BorderLayout.NORTH);
                String swing = option.getAttribute("gui");
                JPanel optionValuesPanel = new JPanel(new FlowLayout());
                optionValuesPanel.setBorder(etchedBorder);
                ButtonGroup bg = null;
                if (swing.equalsIgnoreCase("CHECKBOXES")) {
                    // CHECKBOXES allows selection of 1 or more values from a
                    // predefined set (Analyze, for example)
                    values = option.getElementsByTagName("Value");
                    for (int j = 0; j < values.getLength(); j++) {
                        value = (Element) values.item(j);
                        JCheckBox jcb = new JCheckBox(value.getAttribute("name"));
                        jcb.addMouseListener(this);
                        jcb.setName(flag);
                        if (defaultOption != null
                                && jcb.getActionCommand().equalsIgnoreCase(
                                        defaultOption)) {
                            jcb.setSelected(true);
                        }
                        optionValuesPanel.add(jcb);
                    }
                } else if (swing.equalsIgnoreCase("TEXTFIELD")) {
                    // TEXTFIELD takes an arbitrary String as input
                    JTextField jtf = new JTextField(20);
                    jtf.addMouseListener(this);
                    jtf.setName(flag);
                    if (defaultOption != null && defaultOption.equals("ATOMS")) {
                        FFXSystem sys = mainPanel.getHierarchy().getActive();
                        if (sys != null) {
                            jtf.setText("" + sys.getAtomList().size());
                        }
                    } else if (defaultOption != null) {
                        jtf.setText(defaultOption);
                    }
                    optionValuesPanel.add(jtf);
                } else if (swing.equalsIgnoreCase("RADIOBUTTONS")) {
                    // RADIOBUTTONS allows one choice from a set of predifined
                    // values
                    bg = new ButtonGroup();
                    values = option.getElementsByTagName("Value");
                    for (int j = 0; j < values.getLength(); j++) {
                        value = (Element) values.item(j);
                        JRadioButton jrb = new JRadioButton(value.getAttribute("name"));
                        jrb.addMouseListener(this);
                        jrb.setName(flag);
                        bg.add(jrb);
                        if (defaultOption != null
                                && jrb.getActionCommand().equalsIgnoreCase(
                                        defaultOption)) {
                            jrb.setSelected(true);
                        }
                        optionValuesPanel.add(jrb);
                    }
                } else if (swing.equalsIgnoreCase("PROTEIN")) {
                    // Protein allows selection of amino acids for the protein
                    // builder
                    optionValuesPanel.setLayout(new BoxLayout(
                            optionValuesPanel, BoxLayout.Y_AXIS));
                    optionValuesPanel.add(Box.createRigidArea(new Dimension(0,
                            5)));
                    optionValuesPanel.add(getAminoAcidPanel());
                    optionValuesPanel.add(Box.createRigidArea(new Dimension(0,
                            5)));
                    acidComboBox.removeAllItems();
                    JButton add = new JButton("Edit");
                    add.setActionCommand("PROTEIN");
                    add.addActionListener(this);
                    add.setAlignmentX(Button.CENTER_ALIGNMENT);
                    JPanel comboPanel = new JPanel(new FlowLayout(
                            FlowLayout.CENTER));
                    comboPanel.add(acidTextField);
                    comboPanel.add(add);
                    optionValuesPanel.add(comboPanel);
                    optionValuesPanel.add(Box.createRigidArea(new Dimension(0,
                            5)));
                    JButton remove = new JButton("Remove");
                    add.setMinimumSize(remove.getPreferredSize());
                    add.setPreferredSize(remove.getPreferredSize());
                    remove.setActionCommand("PROTEIN");
                    remove.addActionListener(this);
                    remove.setAlignmentX(Button.CENTER_ALIGNMENT);
                    comboPanel = new JPanel(new FlowLayout(FlowLayout.CENTER));
                    comboPanel.add(acidComboBox);
                    comboPanel.add(remove);
                    optionValuesPanel.add(comboPanel);
                    optionValuesPanel.add(Box.createRigidArea(new Dimension(0,
                            5)));
                    optionValuesPanel.add(acidScrollPane);
                    optionValuesPanel.add(Box.createRigidArea(new Dimension(0,
                            5)));
                    JButton reset = new JButton("Reset");
                    reset.setActionCommand("PROTEIN");
                    reset.addActionListener(this);
                    reset.setAlignmentX(Button.CENTER_ALIGNMENT);
                    optionValuesPanel.add(reset);
                    optionValuesPanel.add(Box.createRigidArea(new Dimension(0,
                            5)));
                    acidTextArea.setText("");
                    acidTextField.setText("");
                } else if (swing.equalsIgnoreCase("NUCLEIC")) {
                    // Nucleic allows selection of nucleic acids for the nucleic
                    // acid builder
                    optionValuesPanel.setLayout(new BoxLayout(
                            optionValuesPanel, BoxLayout.Y_AXIS));
                    optionValuesPanel.add(Box.createRigidArea(new Dimension(0,
                            5)));
                    optionValuesPanel.add(getNucleicAcidPanel());
                    optionValuesPanel.add(Box.createRigidArea(new Dimension(0,
                            5)));
                    acidComboBox.removeAllItems();
                    JButton add = new JButton("Edit");
                    add.setActionCommand("NUCLEIC");
                    add.addActionListener(this);
                    add.setAlignmentX(Button.CENTER_ALIGNMENT);
                    JPanel comboPanel = new JPanel(new FlowLayout(
                            FlowLayout.CENTER));
                    comboPanel.add(acidTextField);
                    comboPanel.add(add);
                    optionValuesPanel.add(comboPanel);
                    optionValuesPanel.add(Box.createRigidArea(new Dimension(0,
                            5)));
                    JButton remove = new JButton("Remove");
                    add.setMinimumSize(remove.getPreferredSize());
                    add.setPreferredSize(remove.getPreferredSize());
                    remove.setActionCommand("NUCLEIC");
                    remove.addActionListener(this);
                    remove.setAlignmentX(Button.CENTER_ALIGNMENT);
                    comboPanel = new JPanel(new FlowLayout(FlowLayout.CENTER));
                    comboPanel.add(acidComboBox);
                    comboPanel.add(remove);
                    optionValuesPanel.add(comboPanel);
                    optionValuesPanel.add(Box.createRigidArea(new Dimension(0,
                            5)));
                    optionValuesPanel.add(acidScrollPane);
                    optionValuesPanel.add(Box.createRigidArea(new Dimension(0,
                            5)));
                    JButton button = new JButton("Reset");
                    button.setActionCommand("NUCLEIC");
                    button.addActionListener(this);
                    button.setAlignmentX(Button.CENTER_ALIGNMENT);
                    optionValuesPanel.add(button);
                    optionValuesPanel.add(Box.createRigidArea(new Dimension(0,
                            5)));
                    acidTextArea.setText("");
                    acidTextField.setText("");
                } else if (swing.equalsIgnoreCase("SYSTEMS")) {
                    // SYSTEMS allows selection of an open system
                    JComboBox<FFXSystem> jcb = new JComboBox<>(
                            mainPanel.getHierarchy().getNonActiveSystems());
                    jcb.setSize(jcb.getMaximumSize());
                    jcb.addActionListener(this);
                    optionValuesPanel.add(jcb);
                }
                // Set up a Conditional for this Option
                if (conditional != null) {
                    isEnabled = false;
                    String conditionalName = conditional.getAttribute("name");
                    String conditionalValues = conditional.getAttribute("value");
                    String cDescription = conditional.getAttribute("description");
                    String cpostProcess = conditional.getAttribute("postProcess");
                    if (conditionalName.toUpperCase().startsWith("KEYWORD")) {
                        optionPanel.setName(conditionalName);
                        String keywords[] = conditionalValues.split(" +");
                        if (activeSystem != null) {
                            Hashtable<String, Keyword> systemKeywords = activeSystem.getKeywords();
                            for (String key : keywords) {
                                if (systemKeywords.containsKey(key.toUpperCase())) {
                                    isEnabled = true;
                                }
                            }
                        }
                    } else if (conditionalName.toUpperCase().startsWith("VALUE")) {
                        isEnabled = true;
                        // Add listeners to the values of this option so
                        // the conditional options can be disabled/enabled.
                        for (int j = 0; j < optionValuesPanel.getComponentCount(); j++) {
                            JToggleButton jtb = (JToggleButton) optionValuesPanel.getComponent(j);
                            jtb.addActionListener(this);
                            jtb.setActionCommand("Conditional");
                        }
                        JPanel condpanel = new JPanel();
                        condpanel.setBorder(etchedBorder);
                        JLabel condlabel = new JLabel(cDescription);
                        condlabel.setEnabled(false);
                        condlabel.setName(conditionalValues);
                        JTextField condtext = new JTextField(10);
                        condlabel.setLabelFor(condtext);
                        if (cpostProcess != null) {
                            condtext.setName(cpostProcess);
                        }
                        condtext.setEnabled(false);
                        condpanel.add(condlabel);
                        condpanel.add(condtext);
                        conditionals.add(condlabel);
                        optionPanel.add(condpanel, BorderLayout.SOUTH);
                    } else if (conditionalName.toUpperCase().startsWith(
                            "REFLECTION")) {
                        String[] condModifiers;
                        if (conditionalValues.equalsIgnoreCase("AltLoc")) {
                            condModifiers = activeSystem.getAltLocations();
                            if (condModifiers != null
                                    && condModifiers.length > 1) {
                                isEnabled = true;
                                bg = new ButtonGroup();
                                for (int j = 0; j < condModifiers.length; j++) {
                                    JRadioButton jrbmi = new JRadioButton(
                                            condModifiers[j]);
                                    jrbmi.addMouseListener(this);
                                    bg.add(jrbmi);
                                    optionValuesPanel.add(jrbmi);
                                    if (j == 0) {
                                        jrbmi.setSelected(true);
                                    }
                                }
                            }
                        } else if (conditionalValues.equalsIgnoreCase("Chains")) {
                            condModifiers = activeSystem.getChainNames();
                            if (condModifiers != null
                                    && condModifiers.length > 0) {
                                isEnabled = true;
                                for (int j = 0; j < condModifiers.length; j++) {
                                    JRadioButton jrbmi = new JRadioButton(
                                            condModifiers[j]);
                                    jrbmi.addMouseListener(this);
                                    bg.add(jrbmi);
                                    optionValuesPanel.add(jrbmi, j);
                                }
                            }
                        }
                    }
                }
                optionPanel.add(optionValuesPanel, BorderLayout.CENTER);
                optionPanel.setPreferredSize(optionPanel.getPreferredSize());
                optionsTabbedPane.addTab(option.getAttribute("name"),
                        optionPanel);
                optionsTabbedPane.setEnabledAt(
                        optionsTabbedPane.getTabCount() - 1, isEnabled);
            }
        }
        optionsTabbedPane.setPreferredSize(optionsTabbedPane.getPreferredSize());
        commandPanel.setLayout(borderLayout);
        commandPanel.add(optionsTabbedPane, BorderLayout.CENTER);
        commandPanel.validate();
        commandPanel.repaint();
        loadLogSettings();
        statusLabel.setText("  " + createCommandInput());
    }

    private void loadLogSettings() {
        String selected = (String) logSettings.getSelectedItem();
        logSettings.removeActionListener(this);
        logSettings.removeAllItems();
        File currentLog;
        String fileName;
        File logDir;
        File systemDir;
        // This implies no Open System
        if (activeSystem == null) {
            fileName = ((String) currentCommandBox.getSelectedItem()).toLowerCase()
                    + ".log";
            currentLog = new File(MainPanel.getPWD() + File.separator
                    + fileName);
            logDir = currentLog.getParentFile();
            systemDir = MainPanel.getPWD();
        } else {
            currentLog = activeSystem.getLogFile();
            fileName = currentLog.getName();
            logDir = currentLog.getParentFile();
            systemDir = activeSystem.getFile().getAbsoluteFile().getParentFile();
        }
        File tempLog;
        File newLog;
        if (logDir == null || !logDir.equals(systemDir)) {
            tempLog = new File(systemDir.getAbsolutePath() + File.separator
                    + fileName);
            if (!commandFileTypes.contains(FileType.ANY)) {
                activeSystem.setLogFile(tempLog);
            }
        } else {
            tempLog = currentLog;
        }
        // Simple Case - default log file doesn't exist yet
        String createNew = null;
        if (!tempLog.exists()) {
            createNew = "Create " + fileName;
            logSettings.addItem(createNew);
            logSettings.setSelectedItem(createNew);
            logString = getLogString(tempLog);
            logSettings.addActionListener(this);
            if (!commandFileTypes.contains(FileType.ANY)) {
                activeSystem.setLogFile(tempLog);
            }
            return;
        }
        // The default log exists, so we can append to it, overwrite it, or
        // create a new one
        newLog = SystemFilter.version(tempLog);
        tempLog = SystemFilter.previousVersion(newLog);
        fileName = tempLog.getName();
        String append = "Append to " + fileName;
        String overwrite = "Overwrite " + fileName;
        logSettings.addItem(append);
        logSettings.addItem(overwrite);
        if (!newLog.equals(tempLog)) {
            createNew = "Create " + newLog.getName();
            logSettings.addItem(createNew);
        }
        if (selected == null) {
            logSettings.setSelectedIndex(0);
            logString = null;
        } else if (selected.startsWith("Append")) {
            logSettings.setSelectedItem(append);
            logString = getLogString(tempLog);
        } else if (selected.startsWith("Overwrite")) {
            logSettings.setSelectedItem(overwrite);
            logString = getLogString(tempLog);
        } else {
            if (createNew != null) {
                logSettings.setSelectedItem(createNew);
                logString = getLogString(newLog);
            } else {
                logString = getLogString(tempLog);
                logSettings.setSelectedItem(append);
            }
        }
        logSettings.addActionListener(this);
    }

    /**
     * {@inheritDoc}
     *
     * Mouse events are used to trigger status bar updates.
     *
     * @param evt
     */
    @Override
    public void mouseClicked(MouseEvent evt) {
        statusLabel.setText("  " + createCommandInput());
    }

    /**
     * {@inheritDoc}
     *
     * @param evt
     */
    @Override
    public void mouseEntered(MouseEvent evt) {
        mouseClicked(evt);
    }

    /**
     * {@inheritDoc}
     *
     * @param evt
     */
    @Override
    public void mouseExited(MouseEvent evt) {
        mouseClicked(evt);
    }

    /**
     * {@inheritDoc}
     *
     * @param evt
     */
    @Override
    public void mousePressed(MouseEvent evt) {
        mouseClicked(evt);
    }

    /**
     * {@inheritDoc}
     *
     * @param evt
     */
    @Override
    public void mouseReleased(MouseEvent evt) {
        mouseClicked(evt);
    }

    /**
     * If a TINKER END file exists for the active system & command, remove it.
     */
    private void removeEnd() {
        FFXSystem m = mainPanel.getHierarchy().getActive();
        if (m == null) {
            return;
        }
        File f = m.getFile();
        if (f == null) {
            return;
        }
        File end = null;
        try {
            if (f.getName().indexOf(".") > 0) {
                String name = f.getName().substring(0,
                        f.getName().lastIndexOf("."));
                end = new File(f.getParent(), name + ".end");
            } else {
                end = new File(f.getParent(), f.getName() + ".end");
            }
        } catch (Exception e) {
        }
        if (end != null && end.exists()) {
            end.delete();
        }
    }

    /**
     * Save ModelingPanel user preferences
     */
    public void savePrefs() {
        String c = ModelingPanel.class.getName();
        prefs.putBoolean(c + ".description", descriptCheckBox.isSelected());
    }

    /**
     **********************************************************************
     */
    // Initialization code and misc. methods.
    public void selected() {
        loadLogSettings();
        setDivider(descriptCheckBox.isSelected());
        validate();
        repaint();
    }

    /**
     * <p>
     * setCommand</p>
     *
     * @param command a {@link java.lang.String} object.
     * @return a boolean.
     */
    public boolean setCommand(String command) {
        if (command == null) {
            return false;
        }
        command = command.toLowerCase();
        command = command.replaceFirst(command.substring(0, 1), command.toUpperCase().substring(0, 1));
        currentCommandBox.setSelectedItem(command);
        mainPanel.setPanel(MainPanel.MODELING);
        return currentCommandBox.getSelectedItem().equals(command);
    }

    /**
     * Set the description divider location
     *
     * @param showDescription True to show the command description
     */
    private void setDivider(boolean showDescription) {
        descriptCheckBox.setSelected(showDescription);
        if (showDescription) {
            int spDivider = (int) (this.getHeight() * (3.0f / 5.0f));
            splitPane.setDividerLocation(spDivider);
        } else {
            splitPane.setDividerLocation(1.0);
        }
    }

    /**
     * Prompt the user to toggle the existence of a TINKER END file.
     */
    private void setEnd() {
        if (ffxThread != null && ffxThread.isAlive()) {
            ModelingShell modelingShell = mainPanel.getModelingShell();
            modelingShell.doInterrupt();
        }
    }

    /**
     **********************************************************************
     *
     * @param mode a {@link java.lang.String} object.
     */
    // Logging Configuation
    public void setLogMode(String mode) {
        mode = mode.toUpperCase();
        for (int i = 0; i < logSettings.getItemCount(); i++) {
            String logType = (String) logSettings.getItemAt(i);
            if (logType.toUpperCase().startsWith(mode)) {
                logSettings.setSelectedIndex(i);
                break;
            }
        }
        loadLogSettings();
    }

    /**
     * <p>
     * deleteLogs</p>
     */
    public void deleteLogs() {
        if (activeSystem != null) {
            File file = activeSystem.getFile();
            String name = file.getName();
            String dir = file.getParent();
            int i = JOptionPane.showConfirmDialog(this, "Delete all logs for "
                    + name + " from " + dir + " ?", "Delete Logs",
                    JOptionPane.YES_NO_OPTION, JOptionPane.WARNING_MESSAGE);
            if (i == JOptionPane.YES_OPTION) {
                try {
                    File files[] = activeSystem.getFile().getParentFile().listFiles();
                    for (File f : files) {
                        name = FilenameUtils.getBaseName(f.getAbsolutePath());
                        if (FilenameUtils.wildcardMatch(f.getName(), name
                                + ".log*")) {
                            f.delete();
                        }
                    }
                } catch (Exception e) {
                }
                activeSystem.setLogFile(null);
                loadLogSettings();
            }
        }
    }

    /**
     * <p>
     * toString</p>
     *
     * @return a {@link java.lang.String} object.
     */
    @Override
    public String toString() {
        return "Modeling Panel";
    }
}
