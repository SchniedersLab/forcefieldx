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

import java.awt.Container;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.logging.Logger;

import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.UIManager;

/**
 * The KeyFileEditor class is a wrapper for the KeywordPanel to create a stand
 * alone Key File Editor (it needs updating).
 *
 * @author schnied
 * @version $Id: $
 */
public final class KeyFileEditor extends JFrame {

    private static final long serialVersionUID = 1L;
    private static final Logger logger = Logger.getLogger(KeyFileEditor.class.getName());

    /**
     * <p>main</p>
     *
     * @param args an array of {@link java.lang.String} objects.
     */
    public static void main(String[] args) {
        KeyFileEditor editor = new KeyFileEditor();
        editor.setVisible(true);
    }
    KeywordPanel tkp;

    /**
     * <p>Constructor for KeyFileEditor.</p>
     */
    public KeyFileEditor() {
        try {
            UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
        } catch (Exception e) {
            logger.warning("Can't set look and feel: " + e);
        }
        tkp = new KeywordPanel(null);
        Container contentPane = getContentPane();
        contentPane.add(tkp);
        setTitle("Key File Editor");
        setSize(800, 800);
        addWindowListener(new WindowAdapter() {

            public void windowClosing(WindowEvent e) {
                close();
                System.exit(0);
            }
        });
    }

    private void close() {
        if (tkp.isFileOpen()) {
            int option = JOptionPane.showConfirmDialog(this,
                    "Save Changes First", "Closing Key Editor",
                    JOptionPane.YES_NO_OPTION, JOptionPane.INFORMATION_MESSAGE);
            if (option == JOptionPane.YES_OPTION) {
                tkp.keySave(null);
            }
        }
    }
}
