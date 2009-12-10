
import javax.swing.*;

JMenuBar menu = cmd.getMainMenu();

JMenu file = menu.getComponent(0);
file.getMenuComponent(0).getAction().putValue(Action.SMALL_ICON, new ImageIcon("/Data/java/icons/folder_page.png"));
file.getMenuComponent(1).getAction().putValue(Action.SMALL_ICON, new ImageIcon("/Data/java/icons/disk.png"));
file.getMenuComponent(2).getAction().putValue(Action.SMALL_ICON, new ImageIcon("/Data/java/icons/cancel.png"));
file.getMenuComponent(5).getAction().putValue(Action.SMALL_ICON, new ImageIcon("/Data/java/icons/drive_web.png"));
// file.remove(10);

JMenu select =  menu.getComponent(1);
select.getMenuComponent(0).getAction().putValue(Action.SMALL_ICON, new ImageIcon("/Data/java/icons/add.png"));
select.getMenuComponent(2).getAction().putValue(Action.SMALL_ICON, new ImageIcon("/Data/java/icons/arrow_merge.png"));
select.getMenuComponent(4).getAction().putValue(Action.SMALL_ICON, new ImageIcon("/Data/java/icons/asterisk_yellow.png"));

JMenu options =  menu.getComponent(4);
options.getMenuComponent(7).getAction().putValue(Action.SMALL_ICON, new ImageIcon("/Data/java/icons/arrow_refresh.png"));
options.getMenuComponent(14).getAction().putValue(Action.SMALL_ICON, new ImageIcon("/Data/java/icons/house.png"));

JMenu pick =  menu.getComponent(5);
pick.getMenuComponent(0).getAction().putValue(Action.SMALL_ICON, new ImageIcon("/Data/java/icons/wand.png"));

JMenu traj =  menu.getComponent(6);
traj.getMenuComponent(0).getAction().putValue(Action.SMALL_ICON, new ImageIcon("/Data/java/icons/control_repeat_blue.png"));
traj.getMenuComponent(5).getAction().putValue(Action.SMALL_ICON, new ImageIcon("/Data/java/icons/control_play_blue.png"));
traj.getMenuComponent(6).getAction().putValue(Action.SMALL_ICON, new ImageIcon("/Data/java/icons/control_stop_blue.png"));
traj.getMenuComponent(7).getAction().putValue(Action.SMALL_ICON, new ImageIcon("/Data/java/icons/control_fastforward_blue.png"));
traj.getMenuComponent(8).getAction().putValue(Action.SMALL_ICON, new ImageIcon("/Data/java/icons/control_rewind_blue.png"));
traj.getMenuComponent(9).getAction().putValue(Action.SMALL_ICON, new ImageIcon("/Data/java/icons/control_start_blue.png"));

JMenu export =  menu.getComponent(8);
export.getMenuComponent(0).getAction().putValue(Action.SMALL_ICON, new ImageIcon("/Data/java/icons/camera.png"));

JMenu window =  menu.getComponent(9);
window.getMenuComponent(0).getAction().putValue(Action.SMALL_ICON, new ImageIcon("/Data/java/icons/application_home.png"));
window.getMenuComponent(5).getAction().putValue(Action.SMALL_ICON, new ImageIcon("/Data/java/icons/application_osx_terminal.png"));

JMenu help =  menu.getComponent(10);
help.getMenuComponent(0).getAction().putValue(Action.SMALL_ICON, new ImageIcon("/Data/java/icons/help.png"));