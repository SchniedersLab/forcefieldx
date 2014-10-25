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
package ffx.potential.bonded;

/**
 * Utilities for creating Nucleic Acid residues.
 *
 * @author Michael Schnieders
 */
public class NucleicAcidUtils {

    /**
     * Biotype keys for nucleic acid backbone atom types. These are consistent
     * with parameter files from TINKER v. 6.1 (June 2012).
     */
    public static final int o5Typ[] = {
        1001, 1031, 1062, 1090, 1117, 1146, 1176, 1203, 0, 0, 0, 0,
        1300, 1334, 1363, 1400, 1431, 1465, 1492, 1523, 1559, 1589, 1624
    };
    public static final int c5Typ[] = {
        1002, 1032, 1063, 1091, 1118, 1147, 1177, 1204, 0, 0, 0, 0,
        1301, 1335, 1364, 1401, 1432, 1466, 1493, 1524, 1560, 1590, 1625
    };
    public static final int h51Typ[] = {
        1003, 1033, 1064, 1092, 1119, 1148, 1178, 1205, 0, 0, 0, 0,
        1302, 1336, 1365, 1402, 1433, 1467, 1494, 1525, 1561, 1591, 1626
    };
    public static final int h52Typ[] = {
        1004, 1034, 1065, 1093, 1120, 1149, 1179, 1206, 0, 0, 0, 0,
        1303, 1337, 1366, 1403, 1434, 1468, 1495, 1526, 1562, 1592, 1627
    };
    public static final int c4Typ[] = {
        1005, 1035, 1066, 1094, 1121, 1150, 1180, 1207, 0, 0, 0, 0,
        1304, 1338, 1367, 1404, 1435, 1469, 1496, 1527, 1563, 1593, 1628
    };
    public static final int h4Typ[] = {
        1006, 1036, 1067, 1095, 1122, 1151, 1181, 1208, 0, 0, 0, 0,
        1305, 1339, 1368, 1405, 1436, 1470, 1497, 1528, 1564, 1594, 1629
    };
    public static final int o4Typ[] = {
        1007, 1037, 1068, 1096, 1123, 1152, 1182, 1209, 0, 0, 0, 0,
        1306, 1340, 1369, 1406, 1437, 1471, 1498, 1529, 1565, 1595, 1630
    };
    public static final int c1Typ[] = {
        1008, 1038, 1069, 1097, 1124, 1153, 1183, 1210, 0, 0, 0, 0,
        1307, 1341, 1370, 1407, 1438, 1472, 1499, 1530, 1566, 1596, 1631
    };
    public static final int h1Typ[] = {
        1009, 1039, 1070, 1098, 1125, 1154, 1184, 1211, 0, 0, 0, 0,
        1308, 1342, 1371, 1408, 1439, 1473, 1500, 1531, 1567, 1597, 1632
    };
    public static final int c3Typ[] = {
        1010, 1040, 1071, 1099, 1126, 1155, 1185, 1212, 0, 0, 0, 0,
        1309, 1343, 1372, 1409, 1440, 1474, 1501, 1532, 1568, 1598, 1633
    };
    public static final int h3Typ[] = {
        1011, 1041, 1072, 1100, 1127, 1156, 1186, 1213, 0, 0, 0, 0,
        1310, 1344, 1373, 1410, 1441, 1475, 1502, 1533, 1569, 1599, 1634
    };
    public static final int c2Typ[] = {
        1012, 1042, 1073, 1101, 1128, 1157, 1187, 1214, 0, 0, 0, 0,
        1311, 1345, 1374, 1411, 1442, 1476, 1503, 1534, 1570, 1600, 1635
    };
    public static final int h21Typ[] = {
        1013, 1043, 1074, 1102, 1129, 1158, 1188, 1215, 0, 0, 0, 0,
        1312, 1346, 1375, 1412, 1443, 1477, 1504, 1535, 1571, 1601, 1636
    };
    public static final int o2Typ[] = {
        1014, 1044, 1075, 1103, 0, 0, 0, 0, 0, 0, 0, 0,
        1313, 1347, 1376, 1413, 1444, 1478, 1505, 1536, 1572, 1602, 1637
    };
    public static final int h22Typ[] = {
        1015, 1045, 1076, 1104, 1130, 1159, 1189, 1216, 0, 0, 0, 0,
        1314, 1348, 1377, 0, 0, 1479, 1506, 1537, 1573, 1603, 1638
    };
    public static final int o3Typ[] = {
        1016, 1046, 1077, 1105, 1131, 1160, 1190, 1217, 0, 0, 0, 0,
        1315, 1349, 1378, 1414, 1445, 1480, 1507, 1538, 1574, 1604, 1639
    };
    public static final int pTyp[] = {
        1230, 1230, 1230, 1230, 1242, 1242, 1242, 1242, 0, 0, 0, 0,
        1230, 1230, 1230, 1230, 1230, 1230, 1230, 1230, 1230, 1230, 1230
    };
    public static final int opTyp[] = {
        1231, 1231, 1231, 1231, 1243, 1243, 1243, 1243, 0, 0, 0, 0,
        1231, 1231, 1231, 1231, 1231, 1231, 1231, 1231, 1231, 1231, 1231
    };
    public static final int h5tTyp[] = {
        1233, 1233, 1233, 1233, 1245, 1245, 1245, 1245, 0, 0, 0, 0,
        1233, 1233, 1233, 1233, 1233, 1233, 1233, 1233, 1233, 1233, 1233
    };
    public static final int h3tTyp[] = {
        1238, 1238, 1238, 1238, 1250, 1250, 1250, 1250, 0, 0, 0, 0,
        1238, 1238, 1238, 1238, 1238, 1238, 1238, 1238, 1238, 1238, 1238
    };

}
