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

import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTextArea;
import javax.swing.JToolBar;
import javax.swing.border.Border;
import javax.swing.border.EtchedBorder;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.prefs.Preferences;

import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.configuration2.Configuration;
import org.apache.commons.configuration2.PropertiesConfiguration;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import ffx.potential.parsers.KeyFileFilter;
import ffx.potential.parsers.KeyFilter;
import ffx.ui.commands.DTDResolver;
import ffx.utilities.Keyword;

/**
 * The KeywordPanel class provides a View and Control of TINKER Keyword (*.KEY)
 * files.
 *
 * @author Michael J. Schnieders
 */
public final class KeywordPanel extends JPanel implements ActionListener {

    private static final Logger logger = Logger.getLogger(KeywordPanel.class.getName());
    private static final Preferences preferences = Preferences.userNodeForPackage(KeywordPanel.class);
    /**
     * True if a Key File is open.
     */
    private boolean fileOpen = false;
    /**
     * Currently open Keywords.
     */
    private Hashtable<String, Keyword> currentKeys;
    /**
     * Currently open Key File.
     */
    private File currentKeyFile;
    /**
     * FFXSystem associated with currently open Key File (if any).
     */
    private FFXSystem currentSystem;
    /**
     * The MainPanel has references to many things the KeywordPanel uses.
     */
    private final MainPanel mainPanel;
    /**
     * HashMap for Keyword GUI Components.
     */
    private LinkedHashMap<String, KeywordComponent> keywordHashMap;
    /**
     * HashMap for Keyword Groups.
     */
    private LinkedHashMap<String, String> groupHashMap;
    /**
     * JComboBox with Keyword Groups, plus a few special cases (Active,
     * FlatFile).
     */
    private JComboBox<String> groupComboBox;
    /**
     * The editPanel holds the toolBar (north), splitPane (center) and
     * statusLabel (south).
     */
    private JPanel editPanel;
    private JToolBar toolBar;
    /**
     * The splitPane holds the editScrollPane (top) and descriptScrollPane
     * (bottom).
     */
    private JSplitPane splitPane;
    private final JLabel statusLabel = new JLabel("  ");
    /**
     * The editScrollPane holds the gridPanel, where KeywordComponents actually
     * live.
     */
    private JScrollPane editScrollPane;
    /**
     * The descriptScrollPane holds the descriptTextArea for Keyword
     * Descriptions.
     */
    private JScrollPane descriptScrollPane;
    /**
     * The descriptTextArea actually holds Keyword Descriptions.
     */
    private JTextArea descriptTextArea;
    /**
     * Allow the user to show/hide Keyword Descriptions.
     */
    private JCheckBoxMenuItem descriptCheckBox;
    private final GridBagLayout gridBagLayout = new GridBagLayout();
    private final GridBagConstraints gridBagConstraints = new GridBagConstraints();
    /**
     * The gridPanel holds an array of KeywordComponents.
     */
    private final JPanel gridPanel = new JPanel(gridBagLayout);
    /**
     * Lines in Keyword files that are comments, unrecognized keywords, or
     * keywords where editing is not supported are stored in a big
     * StringBuilder.
     */
    private StringBuilder commentStringBuffer = new StringBuilder();
    /**
     * This component shows what the saved Key file will look like (WYSIWYG).
     */
    private JTextArea flatfileTextArea;
    /**
     * A simple label if no Keyword File is open.
     */
    private final JLabel noSystemLabel = new JLabel(
            "Keywords for the active system are edited here. ");
    // A simple label if no Keyword Description is available.
    private final JLabel noKeywordLabel = new JLabel(
            "Keyword desciptions are displayed here.");
    private final JPanel noKeywordPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 5, 5));
    private final FlowLayout flowLayout = new FlowLayout(FlowLayout.LEFT, 5, 5);
    private final BorderLayout borderLayout = new BorderLayout();
    private String[] paramNames = null;
    private LinkedHashMap<String, String> paramHashtable = null;

    /**
     * Default Construtor where parent is the tinker Window Frame object.
     *
     * @param f a {@link ffx.ui.MainPanel} object.
     */
    KeywordPanel(MainPanel f) {
        super();
        mainPanel = f;
        initialize();
    }

    /**
     * Handles input from KeywordPanel ToolBar buttons.
     * <p>
     * {@inheritDoc}
     */
    @Override
    public void actionPerformed(ActionEvent evt) {
        String arg = evt.getActionCommand();
        switch (arg) {
            case "Open...":
                keyOpen();
                break;
            case "Close":
                keyClose();
                break;
            case "Save":
                keySave(currentKeyFile);
                break;
            case "Save As...":
                keySaveAs();
                break;
            case "Description":
                JCheckBoxMenuItem box = (JCheckBoxMenuItem) evt.getSource();
                setDivider(box.isSelected());
                break;
        }
        if (evt.getSource() instanceof JComboBox) {
            loadKeywordGroup();
        }
    }

    /**
     * <p>
     * getKeyword</p>
     *
     * @param key a {@link java.lang.String} object.
     * @return a {@link java.lang.String} object.
     */
    public String getKeyword(String key) {
        if (key == null) {
            return null;
        }
        synchronized (this) {
            KeywordComponent keyword = keywordHashMap.get(key.toUpperCase());
            if (keyword == null) {
                return null;
            }
            return keyword.toString();
        }
    }

    /**
     * <p>
     * getKeywordDescription</p>
     *
     * @param key a {@link java.lang.String} object.
     * @return a {@link java.lang.String} object.
     */
    public String getKeywordDescription(String key) {
        if (key == null) {
            return null;
        }
        synchronized (this) {
            KeywordComponent keyword = keywordHashMap.get(key.toUpperCase());
            if (keyword == null) {
                return null;
            }
            return keyword.getKeywordDescription();
        }
    }

    /**
     * <p>
     * getKeywordValue</p>
     *
     * @param key a {@link java.lang.String} object.
     * @return a {@link java.lang.String} object.
     */
    public String getKeywordValue(String key) {
        if (key == null) {
            return null;
        }
        synchronized (this) {
            KeywordComponent keyword = keywordHashMap.get(key.toUpperCase());
            if (keyword == null) {
                return null;
            }
            String value = keyword.toString().trim();
            int firstSpace = value.indexOf(" ");
            if (firstSpace < 1) {
                return value;
            }
            return value.substring(firstSpace, value.length());
        }
    }

    /**
     * <p>
     * Getter for the field <code>paramFiles</code>.</p>
     *
     * @return an array of {@link java.lang.String} objects.
     */
    String[] getParamFiles() {
        return paramNames;
    }

    /**
     * <p>
     * getParamPath</p>
     *
     * @param key a {@link java.lang.String} object.
     * @return a {@link java.lang.String} object.
     */
    public String getParamPath(String key) {
        return paramHashtable.get(key);
    }

    private void initialize() {
        // Load the Keyword Definitions
        loadXML();
        // TextAreas
        flatfileTextArea = new JTextArea();
        flatfileTextArea.setEditable(false);
        flatfileTextArea.setFont(Font.decode("monospaced plain 12"));
        Insets insets = flatfileTextArea.getInsets();
        insets.set(5, 5, 5, 5);
        flatfileTextArea.setMargin(insets);
        // Keyword Edit Panel
        editPanel = new JPanel(flowLayout);
        ClassLoader loader = getClass().getClassLoader();
        ImageIcon icKeyPanel = new ImageIcon(loader.getResource("ffx/ui/icons/page_key.png"));
        noSystemLabel.setIcon(icKeyPanel);
        ImageIcon icon = new ImageIcon(loader.getResource("ffx/ui/icons/information.png"));
        noKeywordLabel.setIcon(icon);
        noKeywordPanel.add(noKeywordLabel);
        editScrollPane = new JScrollPane(editPanel,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        descriptScrollPane = new JScrollPane(descriptTextArea,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        Border eb = BorderFactory.createEtchedBorder(EtchedBorder.RAISED);
        descriptScrollPane.setBorder(eb);
        // Add the Keyword Group Panel and Decription Panel to a JSplitPane
        splitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, editScrollPane,
                descriptScrollPane);
        splitPane.setResizeWeight(1.0);
        splitPane.setOneTouchExpandable(true);
        statusLabel.setBorder(eb);
        // Add the main pieces to the overall KeywordPanel (except the ToolBar)
        setLayout(new BorderLayout());
        add(splitPane, BorderLayout.CENTER);
        add(statusLabel, BorderLayout.SOUTH);
        // Init the GridBagContraints
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 0;
        gridBagConstraints.anchor = GridBagConstraints.WEST;
        gridBagConstraints.gridheight = 1;
        gridBagConstraints.gridwidth = 1;
        gridBagConstraints.fill = GridBagConstraints.HORIZONTAL;
        initToolBar();
        add(toolBar, BorderLayout.NORTH);
        setParamPath();
        loadPrefs();
        loadKeywordGroup();
    }

    /**
     * <p>
     * initToolBar</p>
     */
    private void initToolBar() {
        toolBar = new JToolBar("Keyword Editor");
        toolBar.setLayout(flowLayout);
        ClassLoader loader = getClass().getClassLoader();
        JButton jbopen = new JButton(new ImageIcon(loader.getResource("ffx/ui/icons/folder_page.png")));
        jbopen.setActionCommand("Open...");
        jbopen.setToolTipText("Open a *.KEY File for Editing");
        jbopen.addActionListener(this);
        Insets insets = jbopen.getInsets();
        insets.top = 2;
        insets.bottom = 2;
        insets.left = 2;
        insets.right = 2;
        jbopen.setMargin(insets);
        //toolBar.add(jbopen);
        JButton jbsave = new JButton(new ImageIcon(loader.getResource("ffx/ui/icons/disk.png")));
        jbsave.setActionCommand("Save");
        jbsave.setToolTipText("Save the Active *.KEY File");
        jbsave.addActionListener(this);
        jbsave.setMargin(insets);
        toolBar.add(jbsave);
        JButton jbsaveas = new JButton(new ImageIcon(loader.getResource("ffx/ui/icons/page_save.png")));
        jbsaveas.setActionCommand("Save As...");
        jbsaveas.setToolTipText("Save the Active *.KEY File Under a New Name");
        jbsaveas.addActionListener(this);
        jbsaveas.setMargin(insets);
        toolBar.add(jbsaveas);
        JButton jbclose = new JButton(new ImageIcon(loader.getResource("ffx/ui/icons/cancel.png")));
        jbclose.setActionCommand("Close");
        jbclose.setToolTipText("Close the Active *.KEY File");
        jbclose.addActionListener(this);
        jbclose.setMargin(insets);
        //toolBar.add(jbclose);
        toolBar.addSeparator();
        groupComboBox.setMaximumSize(groupComboBox.getPreferredSize());
        groupComboBox.addActionListener(this);
        groupComboBox.setEditable(false);
        toolBar.add(groupComboBox);
        toolBar.addSeparator();
        ImageIcon icinfo = new ImageIcon(loader.getResource("ffx/ui/icons/information.png"));
        descriptCheckBox = new JCheckBoxMenuItem(icinfo);
        descriptCheckBox.setActionCommand("Description");
        descriptCheckBox.setToolTipText("Show/Hide Keyword Descriptions");
        descriptCheckBox.addActionListener(this);
        descriptCheckBox.setMargin(insets);
        toolBar.add(descriptCheckBox);
        toolBar.add(new JLabel(""));
        toolBar.setBorderPainted(false);
        toolBar.setFloatable(false);
        toolBar.setRollover(true);
        toolBar.setOrientation(JToolBar.HORIZONTAL);
    }

    /**
     * <p>
     * isFileOpen</p>
     *
     * @return a boolean.
     */
    boolean isFileOpen() {
        return fileOpen;
    }

    /**
     * <p>
     * isKeyword</p>
     *
     * @param key a {@link java.lang.String} object.
     * @return a boolean.
     */
    public boolean isKeyword(String key) {
        synchronized (this) {
            KeywordComponent keyword = keywordHashMap.get(key.toUpperCase());
            return keyword != null;
        }
    }

    /**
     * <p>
     * keyClear</p>
     */
    private void keyClear() {
        synchronized (this) {
            // Clear each KeywordComponent
            for (KeywordComponent kw : keywordHashMap.values()) {
                kw.clearKeywordComponent();
            }
            KeywordComponent.setKeywordModified(false);
            // Reset internal variables
            fileOpen = false;
            currentKeys = null;
            currentSystem = null;
            currentKeyFile = null;
            // Reset the View
            commentStringBuffer = new StringBuilder();
            statusLabel.setText("  ");
            loadKeywordGroup();
        }
    }

    /**
     * <p>
     * keyClose</p>
     *
     * @return a boolean.
     */
    private boolean keyClose() {
        if (KeywordComponent.isKeywordModified()) {
            return false;
        }
        keyClear();
        return true;
    }

    /**
     * Give the user a File Dialog Box so they can select a key file.
     */
    private void keyOpen() {
        if (fileOpen && KeywordComponent.isKeywordModified()) {
            int option = JOptionPane.showConfirmDialog(this,
                    "Save Changes First", "Opening New File",
                    JOptionPane.YES_NO_CANCEL_OPTION,
                    JOptionPane.INFORMATION_MESSAGE);
            if (option == JOptionPane.CANCEL_OPTION) {
                return;
            } else if (option == JOptionPane.YES_OPTION) {
                keySave(currentKeyFile);
            }
            keyClear();
        }
        JFileChooser d = MainPanel.resetFileChooser();
        if (currentSystem != null) {
            File cwd = currentSystem.getFile();
            if (cwd != null && cwd.getParentFile() != null) {
                d.setCurrentDirectory(cwd.getParentFile());
            }
        }
        d.setAcceptAllFileFilterUsed(false);
        d.setFileFilter(MainPanel.keyFileFilter);
        d.setDialogTitle("Open KEY File");
        int result = d.showOpenDialog(this);
        if (result == JFileChooser.APPROVE_OPTION) {
            File newKeyFile = d.getSelectedFile();
            if (newKeyFile != null && newKeyFile.exists()
                    && newKeyFile.canRead()) {
                keyOpen(newKeyFile);
            }
        }
    }

    /**
     * <p>
     * keyOpen</p>
     *
     * @param newKeyFile a {@link java.io.File} object.
     * @return a boolean.
     */
    private boolean keyOpen(File newKeyFile) {
        if (newKeyFile != null && newKeyFile.exists() && newKeyFile.canRead()) {
            Hashtable<String, Keyword> newKeys = KeyFilter.open(newKeyFile);
            if (newKeys != null && newKeys.size() > 0) {
                if (currentSystem != null) {
                    currentSystem.setKeyFile(currentKeyFile);
                    currentSystem.setKeywords(currentKeys);
                }
                loadActive(currentSystem, newKeys, newKeyFile);
                return true;
            }
        }
        return false;
    }

    /**
     * <p>
     * keySave</p>
     *
     * @param f a {@link java.io.File} object.
     */
    void keySave(File f) {
        if (f != null) {
            currentKeyFile = f;
        }
        if (!fileOpen || currentKeyFile == null) {
            return;
        }
        storeActive();
        saveKeywords(currentKeyFile, keywordHashMap, commentStringBuffer);
    }

    /**
     * <p>
     * keySaveAs</p>
     */
    private void keySaveAs() {
        if (!fileOpen) {
            return;
        }
        JFileChooser d = MainPanel.resetFileChooser();
        d.setDialogTitle("Save KEY File");
        d.setAcceptAllFileFilterUsed(false);
        if (currentKeyFile != null) {
            d.setCurrentDirectory(currentKeyFile.getParentFile());
            d.setSelectedFile(currentKeyFile);
        }
        d.setFileFilter(new KeyFileFilter());
        int result = d.showSaveDialog(this);
        if (result == JFileChooser.APPROVE_OPTION) {
            currentKeyFile = d.getSelectedFile();
            keySave(currentKeyFile);
        }
    }

    /**
     * <p>
     * loadActive</p>
     *
     * @param newSystem a {@link ffx.ui.FFXSystem} object.
     * @return a boolean.
     */
    boolean loadActive(FFXSystem newSystem) {
        synchronized (this) {
            if (newSystem == null) {
                keyClear();
                return false;
            }
            configToKeywords(newSystem);
            File newKeyFile = newSystem.getKeyFile();
            if (newKeyFile != null) {
                //logger.info(String.format("Key File %s.", newKeyFile.getAbsolutePath()));
            }
            Hashtable<String, Keyword> newKeys = newSystem.getKeywords();
            if (newKeyFile == null && newKeys == null) {
                logger.info(String.format("Loaded %s with no keywords.", newSystem.toString()));
                return false;
            }
            //logger.info(String.format("Loading %s with %d keywords.", newSystem.toString(), newKeys.size()));
            return loadActive(newSystem, newKeys, newKeyFile);
        }
    }

    private void configToKeywords(FFXSystem newSystem) {

        CompositeConfiguration properties = newSystem.getProperties();
        Hashtable<String, Keyword> keywordHash = new Hashtable<>();

        // Create the "Comments" property.
        Keyword keyword = new Keyword("COMMENTS");
        keywordHash.put("COMMENTS", keyword);

        // Loop over properties from the keyword file.
        Configuration config = properties.getConfiguration(0);
        if (config instanceof PropertiesConfiguration) {
            PropertiesConfiguration prop = (PropertiesConfiguration) config;
            Iterator<String> keys = prop.getKeys();
            while (keys.hasNext()) {
                String key = keys.next();
                if (keywordHash.contains(key)) {
                    keyword = keywordHash.get(key);
                    keyword.append(prop.getStringArray(key));
                } else {
                    String[] values = prop.getStringArray(key);
                    keyword = new Keyword(key, values);
                    keywordHash.put(key, keyword);
                }
            }
        }

        newSystem.setKeywords(keywordHash);
    }

    /**
     * <p>
     * loadActive</p>
     *
     * @param newSystem  a {@link ffx.ui.FFXSystem} object.
     * @param newKeys    a {@link java.util.Hashtable} object.
     * @param newKeyFile a {@link java.io.File} object.
     * @return a boolean.
     */
    private boolean loadActive(FFXSystem newSystem,
                               Hashtable<String, Keyword> newKeys, File newKeyFile) {

        synchronized (this) {
            // Store changes made to the current system (if any) first.
            if (currentKeys != null && KeywordComponent.isKeywordModified()) {
                if (currentSystem != null && currentSystem != newSystem) {
                    storeActive();
                } else {
                    saveChanges();
                }
            }
            // Clear previous values
            keyClear();
            // No keywords to load
            if (newKeys == null || newKeys.size() == 0) {
                return false;
            }
            currentSystem = newSystem;
            currentKeyFile = newKeyFile;
            currentKeys = newKeys;
            fileOpen = true;
            /*
             * Keys to remove are those entries not recognized by the FFX
             * Keyword Editor These keys are removed from the active System's
             * keyword Hashtable and appended to the generic "Comments"
             * KeywordData instance.
             */
            ArrayList<String> keysToRemove = new ArrayList<>();
            Keyword comments = currentKeys.get("COMMENTS");
            for (Keyword keyword : currentKeys.values()) {
                String label = keyword.getKeyword();
                Vector<String> data = keyword.getEntries();
                if (label.equals("COMMENTS")) {
                    continue;
                }
                KeywordComponent tk = keywordHashMap.get(label.toUpperCase());
                // If the keyword label isn't recognized, preserve the line
                // as a Comment
                if (tk == null) {
                    keysToRemove.add(keyword.getKeyword());
                    if (data.isEmpty()) {
                        comments.append(label);
                    } else if (label.equalsIgnoreCase("MULTIPOLE")) {
                        int count = 5;
                        for (String entry : data) {
                            count++;
                            if (count > 5) {
                                comments.append(label + " " + entry);
                                count = 1;
                            } else {
                                comments.append(entry);
                            }
                        }
                    } else if (label.equalsIgnoreCase("TORTORS")) {
                        int points = 0;
                        int count = 0;
                        for (String entry : data) {
                            count++;
                            if (count > points) {
                                String[] res = entry.split(" +");
                                int xres = Integer.parseInt(res[5]);
                                int yres = Integer.parseInt(res[5]);
                                points = xres * yres;
                                comments.append(label + " " + entry);
                                count = 0;
                            } else {
                                comments.append(entry);
                            }
                        }
                    } else {
                        for (String entry : data) {
                            comments.append(label + " " + entry);
                        }
                    }
                    continue;
                }
                // No data for this keyword (something like "verbose" that
                // doesn't have modifiers)
                if (data.isEmpty()) {
                    tk.loadKeywordEntry(label);
                }
                // One to many data entries (something "atom" or "bond" that
                // can be repeated)
                for (String s : data) {
                    tk.loadKeywordEntry(s);
                }
            }
            // Load up the Comments
            Vector<String> entries = comments.getEntries();
            if (entries != null) {
                for (String s : entries) {
                    if (s.length() > 1) {
                        if (!groupHashMap.containsKey(s.substring(1).toUpperCase().trim())) {
                            commentStringBuffer.append(s).append("\n");
                        }
                    } else {
                        commentStringBuffer.append(s).append("\n");
                    }
                }
            }
            // Keywords that aren't recognized have been turned into Comments
            for (String k : keysToRemove) {
                currentKeys.remove(k);
            }
            if (currentSystem != null) {
                currentSystem.setKeywords(currentKeys);
                currentSystem.setKeyFile(currentKeyFile);
            }
            loadKeywordGroup();
            return true;
        }
    }

    private void loadKeywordGroup() {
        synchronized (this) {
            editPanel.removeAll();
            String selectedGroup = (String) groupComboBox.getSelectedItem();
            if (currentKeys == null) {
                editPanel.setLayout(flowLayout);
                editPanel.add(noSystemLabel);
                int temp = splitPane.getDividerLocation();
                splitPane.setBottomComponent(noKeywordPanel);
                splitPane.setDividerLocation(temp);
            } else if (selectedGroup.equalsIgnoreCase("Flat File View")) {
                editPanel.setLayout(borderLayout);
                publishKeywords();
                editPanel.add(flatfileTextArea, BorderLayout.CENTER);
                int temp = splitPane.getDividerLocation();
                splitPane.setBottomComponent(noKeywordPanel);
                splitPane.setDividerLocation(temp);
            } else if (selectedGroup.equalsIgnoreCase("Active Keywords")) {
                gridPanel.removeAll();
                gridBagConstraints.gridy = 0;
                gridBagConstraints.anchor = GridBagConstraints.WEST;
                gridBagConstraints.gridheight = 1;
                gridBagConstraints.gridwidth = 1;
                gridBagConstraints.fill = GridBagConstraints.HORIZONTAL;
                for (KeywordComponent keywordComponent : keywordHashMap.values()) {
                    if (keywordComponent.isActive()) {
                        JPanel jptemp = keywordComponent.getKeywordGUI();
                        gridBagLayout.setConstraints(jptemp, gridBagConstraints);
                        gridPanel.add(jptemp);
                        gridBagConstraints.gridy++;
                    }
                }
                KeywordComponent.fillPanel(gridPanel, gridBagLayout,
                        gridBagConstraints);
                editPanel.setLayout(flowLayout);
                editPanel.add(gridPanel);
                int temp = splitPane.getDividerLocation();
                splitPane.setBottomComponent(descriptScrollPane);
                splitPane.setDividerLocation(temp);
            } else {
                gridPanel.removeAll();
                gridBagConstraints.gridy = 0;
                gridBagConstraints.gridx = 0;
                gridBagConstraints.anchor = GridBagConstraints.WEST;
                gridBagConstraints.gridheight = 1;
                gridBagConstraints.gridwidth = 1;
                gridBagConstraints.fill = GridBagConstraints.HORIZONTAL;
                for (KeywordComponent keywordComponent : keywordHashMap.values()) {
                    if (keywordComponent.getKeywordGroup().equalsIgnoreCase(
                            selectedGroup)) {
                        JPanel jptemp = keywordComponent.getKeywordGUI();
                        gridBagLayout.setConstraints(jptemp, gridBagConstraints);
                        gridPanel.add(jptemp);
                        gridBagConstraints.gridy++;
                    }
                }
                KeywordComponent.fillPanel(gridPanel, gridBagLayout,
                        gridBagConstraints);
                editPanel.setLayout(flowLayout);
                editPanel.add(gridPanel);
                int temp = splitPane.getDividerLocation();
                splitPane.setBottomComponent(descriptScrollPane);
                splitPane.setDividerLocation(temp);
            }
            if (currentKeyFile != null) {
                statusLabel.setText("  " + currentKeyFile.toString());
            } else {
                statusLabel.setText("  ");
            }
            editScrollPane.validate();
            editScrollPane.repaint();
        }
    }

    /**
     * <p>
     * loadPrefs</p>
     */
    private void loadPrefs() {
        String c = KeywordPanel.class.getName();
        descriptCheckBox.setSelected(!preferences.getBoolean(c + ".description", true));
        descriptCheckBox.doClick();
    }

    /**
     * Load up the Force Field X Keyword Definitions
     */
    private void loadXML() {
        NodeList groups, keywords, values;
        Element group, keyword, value;
        String groupName;
        String keywordName, keywordDescription, keywordGUI;
        groups = null;
        try {
            // Build the Document from the keywords.xml file
            DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
            DocumentBuilder db = dbf.newDocumentBuilder();
            db.setEntityResolver(new DTDResolver());
            URL keyURL = getClass().getClassLoader().getResource(
                    "ffx/ui/commands/keywords.xml");
            Document doc = db.parse(keyURL.openStream());
            Element document = doc.getDocumentElement();
            Element body = (Element) document.getElementsByTagName("body").item(0);
            groups = body.getElementsByTagName("section");
        } catch (ParserConfigurationException | SAXException | IOException e) {
            logger.warning(e.toString());
        }
        keywordHashMap = new LinkedHashMap<>();
        groupHashMap = new LinkedHashMap<>();
        groupHashMap.put("ACTIVE KEYWORDS", "Active Keywords");
        groupHashMap.put("FLAT FILE VIEW", "Flat File View");
        groupComboBox = new JComboBox<>();
        groupComboBox.addItem("Active Keywords");
        groupComboBox.addItem("Flat File View");
        descriptTextArea = new JTextArea();
        descriptTextArea.setLineWrap(true);
        descriptTextArea.setWrapStyleWord(true);
        Insets insets = descriptTextArea.getInsets();
        insets.set(5, 5, 5, 5);
        descriptTextArea.setMargin(insets);
        int length = groups.getLength();
        // Iterate through the Keyword Groups
        for (int i = 0; i < length; i++) {
            group = (Element) groups.item(i);
            groupName = group.getAttribute("name");
            groupComboBox.addItem(groupName);
            keywords = group.getElementsByTagName("subsection");
            int klength = keywords.getLength();
            for (int j = 0; j < klength; j++) {
                keyword = (Element) keywords.item(j);
                keywordName = keyword.getAttribute("name");
                Node text = keyword.getFirstChild();
                keywordDescription = text.getNodeValue().replace('\n', ' ');
                keywordDescription = keywordDescription.replaceAll("  +", " ");
                keywordGUI = keyword.getAttribute("rep");
                KeywordComponent.SwingRepresentation type;
                try {
                    type = KeywordComponent.SwingRepresentation.valueOf(keywordGUI.toUpperCase());
                } catch (Exception e) {
                    type = null;
                    logger.log(Level.WARNING, "{0}: Unknown GUI Component - {1}", new Object[]{keywordName, type});
                    System.exit(-1);
                }
                KeywordComponent key;
                if (type == KeywordComponent.SwingRepresentation.CHECKBOXES
                        || type == KeywordComponent.SwingRepresentation.COMBOBOX) {
                    values = keyword.getElementsByTagName("Value");
                    String labels[] = new String[values.getLength()];
                    for (int k = 0; k < values.getLength(); k++) {
                        value = (Element) values.item(k);
                        labels[k] = value.getAttribute("name");
                    }
                    key = new KeywordComponent(keywordName, groupName, type,
                            keywordDescription, descriptTextArea, labels);
                } else {
                    key = new KeywordComponent(keywordName, groupName, type,
                            keywordDescription, descriptTextArea);
                }
                keywordHashMap.put(keywordName.toUpperCase(), key);
                groupHashMap.put(groupName.toUpperCase(), groupName);
            }
        }
        groupComboBox.setSelectedIndex(0);
    }

    /**
     * <p>
     * publishKeywords</p>
     */
    private void publishKeywords() {
        synchronized (this) {
            flatfileTextArea.setText("");
            boolean writegroup = false;
            String pgroup = null;
            // Write out keywords in groups
            for (KeywordComponent keyword : keywordHashMap.values()) {
                String group = keyword.getKeywordGroup();
                if (pgroup == null || !group.equalsIgnoreCase(pgroup)) {
                    writegroup = true;
                    pgroup = group;
                }
                String line = keyword.toString();
                if (line != null) {
                    if (writegroup) {
                        flatfileTextArea.append("\n");
                        flatfileTextArea.append("# " + group);
                        flatfileTextArea.append("\n");
                        writegroup = false;
                    }
                    flatfileTextArea.append(line);
                    flatfileTextArea.append("\n");
                }
            }
            flatfileTextArea.append("\n");
            String s = commentStringBuffer.toString();
            if (!s.trim().equals("")) {
                flatfileTextArea.append(s.trim());
            }
            flatfileTextArea.append("\n");
        }
    }

    /**
     * <p>
     * saveChanges</p>
     *
     * @return a boolean.
     */
    boolean saveChanges() {
        if (KeywordComponent.isKeywordModified() && currentKeyFile != null) {
            keySave(currentKeyFile);
            return true;
        }
        return false;
    }

    /**
     * <p>
     * saveKeywords</p>
     *
     * @param keyFile        a {@link java.io.File} object.
     * @param keywordHashMap a {@link java.util.LinkedHashMap} object.
     * @param comments       a {@link java.lang.StringBuilder} object.
     * @return a boolean.
     */
    private boolean saveKeywords(File keyFile,
                                 LinkedHashMap<String, KeywordComponent> keywordHashMap,
                                 StringBuilder comments) {
        synchronized (this) {
            FileWriter fw = null;
            BufferedWriter bw = null;
            try {
                fw = new FileWriter(keyFile);
                bw = new BufferedWriter(fw);
                boolean writegroup = false;
                String pgroup = null;
                // Write out keywords in groups
                for (KeywordComponent keyword : keywordHashMap.values()) {
                    String group = keyword.getKeywordGroup();
                    if (pgroup == null || !group.equalsIgnoreCase(pgroup)) {
                        writegroup = true;
                        pgroup = group;
                    }
                    String line = keyword.toString();
                    if (line != null) {
                        if (writegroup) {
                            bw.newLine();
                            bw.write("# " + group);
                            bw.newLine();
                            writegroup = false;
                        }
                        bw.write(line);
                        bw.newLine();
                    }
                }
                bw.newLine();
                String s = comments.toString();
                if (s != null && !s.trim().equals("")) {
                    bw.write(s.trim());
                }
                bw.newLine();
                bw.flush();
                KeywordComponent.setKeywordModified(false);
            } catch (IOException e) {
                logger.warning(e.toString());
                return false;
            } finally {
                try {
                    if (bw != null) {
                        bw.close();
                    }
                    if (fw != null) {
                        fw.close();
                    }
                } catch (Exception e) {
                    logger.warning(e.toString());
                }
            }
            return true;
        }
    }

    /**
     * <p>
     * savePrefs</p>
     */
    void savePrefs() {
        String c = KeywordPanel.class.getName();
        preferences.putInt(c + ".divider", splitPane.getDividerLocation());
        preferences.putBoolean(c + ".description", descriptCheckBox.isSelected());
    }

    /**
     * <p>
     * selected</p>
     */
    public void selected() {
        setDivider(descriptCheckBox.isSelected());
        validate();
        repaint();
    }

    /**
     * <p>
     * setDivider</p>
     *
     * @param b a boolean.
     */
    private void setDivider(boolean b) {
        if (b) {
            descriptCheckBox.setSelected(b);
            int spDivider = (int) (this.getHeight() * (3.0f / 5.0f));
            splitPane.setDividerLocation(spDivider);
        } else {
            splitPane.setDividerLocation(1.0);
        }
    }

    /**
     * Make the passed Keyword Group active in the editor.
     *
     * @param keygroup String
     */
    public void setKeywordGroup(String keygroup) {
        synchronized (this) {
            keygroup = groupHashMap.get(keygroup.toUpperCase());
            if (keygroup == null) {
                return;
            }
            if (!groupComboBox.getSelectedItem().equals(keygroup)) {
                groupComboBox.setSelectedItem(keygroup);
                loadKeywordGroup();
            }
        }
    }

    /**
     * Load a value into a KeywordComponent. Value is equivalent to one line in
     * a TINKER key file, except without the keyword at the beginning. Value
     * should be null to just indicate the Keyword is present (active). If this
     * Keyword can apprear many times, value will be appended to the list.
     *
     * @param key   String
     * @param value String
     */
    public void setKeywordValue(String key, String value) {
        synchronized (this) {
            KeywordComponent keyword = keywordHashMap.get(key.toUpperCase());
            if (keyword == null) {
                return;
            }
            keyword.loadKeywordEntry(value);
            String keygroup = keyword.getKeywordGroup();
            if (!groupComboBox.getSelectedItem().equals(keygroup)) {
                groupComboBox.setSelectedItem(keygroup);
                loadKeywordGroup();
            }
            mainPanel.setPanel(MainPanel.KEYWORDS);
        }
    }

    /**
     * <p>
     * setParamPath</p>
     */
    private void setParamPath() {
        // Location of TINKER parameter files
        File paramDir = new File(MainPanel.ffxDir.getAbsolutePath() + File.separator
                + "params");
        if (paramDir.exists()) {
            File[] paramFiles = paramDir.listFiles();
            paramHashtable = new LinkedHashMap<>();
            for (File f : paramFiles) {
                if (f.exists() && f.canRead()
                        && MainPanel.forceFieldFileFilter.accept(f)) {
                    paramHashtable.put(f.getName(), f.getAbsolutePath());
                }
            }
            int num = paramHashtable.size();
            paramNames = new String[num + 1];
            int i = 1;
            for (String name : paramHashtable.keySet()) {
                paramNames[i] = name;
                i++;
            }
            paramNames[0] = "AAA";
            java.util.Arrays.sort(paramNames);
            paramNames[0] = "Use an existing TINKER Key file".intern();
        }
    }

    /**
     * Store the KeywordPanel's current keyword content into sys.
     *
     * @param sys FFXSystem
     */
    public void store(FFXSystem sys) {
        synchronized (this) {
            FFXSystem back = currentSystem;
            currentSystem = sys;
            storeActive();
            currentSystem = back;
        }
    }

    /**
     * Store the KeywordPanel's current keyword content in the activeSystem.
     */
    private void storeActive() {
        synchronized (this) {
            if (currentSystem == null) {
                return;
            }
            // No changes
            if (!KeywordComponent.isKeywordModified()) {
                return;
            }
            Hashtable<String, Keyword> currentKeys = currentSystem.getKeywords();
            Hashtable<String, Keyword> newKeys = new Hashtable<String, Keyword>();
            for (KeywordComponent kc : keywordHashMap.values()) {
                if (kc.isActive()) {
                    Keyword keywordData = null;
                    if (currentKeys != null) {
                        keywordData = currentKeys.get(kc.getKeyword());
                    }
                    if (keywordData == null) {
                        keywordData = new Keyword(kc.getKeyword());
                    } else {
                        keywordData.clear();
                    }
                    kc.getKeywordData(keywordData);
                    newKeys.put(kc.getKeyword(), keywordData);
                }
            }
            Keyword comments = new Keyword("COMMENTS", commentStringBuffer.toString());
            newKeys.put("COMMENTS", comments);
            currentSystem.setKeywords(newKeys);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        return "Keyword Editor";
    }
}
