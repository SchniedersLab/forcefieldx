/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.ui;

import java.awt.Dimension;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Hashtable;
import java.util.logging.Logger;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JSlider;
import javax.swing.border.Border;
import javax.swing.border.EtchedBorder;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import ffx.potential.bonded.MSRoot;
import ffx.potential.bonded.RendererCache;

/**
 * The GraphicsPrefs class allows users to select graphics preferences.
 *
 * @author schnied
 * @version $Id: $
 */
public class GraphicsPrefs extends JDialog implements ActionListener {

    private static final Logger logger = Logger.getLogger(GraphicsPrefs.class.getName());
    private static final long serialVersionUID = 1L;
    private MSRoot root;
    private GridBagConstraints constraints;
    private boolean change = false;

    /**
     * Contructor
     *
     * @param frame
     *            Parent frame
     * @param r
     *            Data structure root
     */
    public GraphicsPrefs(Frame frame, MSRoot r) {
        super(frame, "", true);
        root = r;
        setTitle("Graphics Preferences");
        setSize(400, 200);
        setResizable(false);
        getContentPane().setLayout(new GridBagLayout());
        constraints = new GridBagConstraints();
        constraints.weighty = 100;
        constraints.gridheight = 1;
        constraints.gridwidth = 2;
        constraints.ipadx = 5;
        constraints.ipady = 5;
        constraints.gridx = 0;
        constraints.gridy = 0;
        constraints.insets = new Insets(5, 5, 5, 5);
        // Slider for radius
        JSlider radius = new JSlider(0, 200, 1);
        radius.setMajorTickSpacing(50);
        radius.setMinorTickSpacing(10);
        radius.setPaintLabels(true);
        radius.setPaintTicks(true);
        radius.setValue((int) (RendererCache.radius * 100));
        Hashtable<Integer, JLabel> labelTable = new Hashtable<Integer, JLabel>();
        labelTable.put(new Integer(1), new JLabel("0%"));
        labelTable.put(new Integer(50), new JLabel("50%"));
        labelTable.put(new Integer(100), new JLabel("100%"));
        labelTable.put(new Integer(150), new JLabel("150%"));
        labelTable.put(new Integer(200), new JLabel("200%"));
        radius.setLabelTable(labelTable);
        addSlider(radius, " Radius", 1);
        // Slider for bondwidth
        JSlider bondwidth = new JSlider(1, 5, 5);
        bondwidth.setMajorTickSpacing(1);
        bondwidth.setMinorTickSpacing(1);
        bondwidth.setPaintLabels(true);
        bondwidth.setPaintTicks(true);
        labelTable = new Hashtable<Integer, JLabel>();
        labelTable.put(new Integer(1), new JLabel("1"));
        labelTable.put(new Integer(3), new JLabel("3"));
        labelTable.put(new Integer(5), new JLabel("5"));
        bondwidth.setLabelTable(labelTable);
        bondwidth.setValue(RendererCache.bondwidth);
        addSlider(bondwidth, " Wireframe Thickness", 3);
        // Slider for detail
        JSlider detail = new JSlider(0, 10, 3);
        detail.setMajorTickSpacing(1);
        detail.setMinorTickSpacing(1);
        detail.setPaintLabels(true);
        detail.setPaintTicks(true);
        detail.setValue(RendererCache.detail);
        labelTable = new Hashtable<Integer, JLabel>();
        labelTable.put(new Integer(0), new JLabel("Performance"));
        labelTable.put(new Integer(10), new JLabel("Quality"));
        detail.setLabelTable(labelTable);
        addSlider(detail, " Detail", 2);
        constraints.gridwidth = 1;
        JButton jb = new JButton("Apply");
        jb.addActionListener(this);
        getContentPane().add(jb, constraints);
        JButton jbclose = new JButton("Close");
        jbclose.addActionListener(this);
        constraints.gridx++;
        getContentPane().add(jbclose, constraints);
        pack();
        Dimension dim = getToolkit().getScreenSize();
        Dimension ddim = getSize();
        setLocation((dim.width - ddim.width) / 2,
                (dim.height - ddim.height) / 2);
    }

    /** {@inheritDoc} */
    public void actionPerformed(ActionEvent e) {
        if (e.getActionCommand().equalsIgnoreCase("Apply")) {
            if (change != false) {
                root.setView(RendererCache.ViewModel.DETAIL, null);
            }
            change = false;
        } else if (e.getActionCommand().equalsIgnoreCase("Close")) {
            if (change != false) {
                root.setView(RendererCache.ViewModel.DETAIL, null);
            }
            change = false;
            dispose();
        }
    }

    /**
     * <p>addSlider</p>
     *
     * @param s a {@link javax.swing.JSlider} object.
     * @param description a {@link java.lang.String} object.
     * @param sliderID a int.
     */
    public void addSlider(JSlider s, String description, final int sliderID) {
        Border eb = BorderFactory.createEtchedBorder(EtchedBorder.LOWERED);
        s.setSnapToTicks(true);
        s.setBorder(new TitledBorder(eb, description));
        s.addChangeListener(new ChangeListener() {

            public void stateChanged(ChangeEvent e) {
                JSlider source = (JSlider) e.getSource();
                if (source.getValueIsAdjusting()) {
                    return;
                }
                int value = source.getValue();
                switch (sliderID) {
                    case 1:
                        if (value < 1) {
                            return;
                        }
                        double temp = value / 100.0d;
                        if (temp != RendererCache.radius) {
                            change = true;
                            RendererCache.radius = temp;
                        }
                        break;
                    case 2:
                        if (RendererCache.detail != value) {
                            change = true;
                            RendererCache.detail = value;
                        }
                        break;
                    case 3:
                        if (RendererCache.bondwidth != value) {
                            change = true;
                            RendererCache.bondwidth = value;
                        }
                        break;
                    default:
                        logger.info("Unknown Slider");
                }
            }
        });
        // add three components into the next row
        constraints.gridx = 0;
        constraints.anchor = GridBagConstraints.WEST;
        constraints.fill = GridBagConstraints.NONE;
        getContentPane().add(s, constraints);
        constraints.gridy++;
    }
}
