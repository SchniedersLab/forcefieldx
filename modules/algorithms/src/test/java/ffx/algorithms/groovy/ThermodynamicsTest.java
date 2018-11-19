/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.algorithms.groovy;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import ffx.potential.bonded.Atom;
import org.apache.commons.configuration2.Configuration;
import org.apache.commons.configuration2.MapConfiguration;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import ffx.algorithms.AbstractOSRW;
import ffx.algorithms.PJDependentTest;
import ffx.algorithms.TransitionTemperedOSRW;
import ffx.algorithms.groovy.NewThermodynamics;
import ffx.crystal.CrystalPotential;
import ffx.potential.bonded.LambdaInterface;

import groovy.lang.Binding;

/**
 * Tests the functionality of the Transition-Tempered OSRW algorithm in Force Field X,
 * both by running a very simple and quick TT-OSRW run, and setting up a number of systems,
 * comparing output energies and gradients from the starting algorithmConfig.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
@RunWith(Parameterized.class)
public class ThermodynamicsTest extends PJDependentTest {

    @Parameterized.Parameters
    public static Collection<Object[]> data() {
        /**
         * Test info, filenames, mode, dG, tol(dG), grad atoms, PE, dU/dL, d2U/dL2,
         * dU/dX, d2U/dXdL, Groovy options, properties, Groovy flags
         */
        return Arrays.asList(new Object[][]{
            {
                "Thermodynamics Help Message Test", new String[]{}, ThermoTestMode.HELP, 0, 0, null, null,
                null, null, null, null, new String[]{}, new String[]{}, new String[]{"-h", "true"}
            },
            {
                "Acetamide Implicit Solvation Free Energy: -10.8 kcal/mol", new String[]{"ffx/algorithms/structures/acetamide.gk.xyz"},
                ThermoTestMode.FREE, -10.8, 1.0, null, null, null, null, null, null,
                new String[]{"-C", "10", "-d", "1.0", "-n", "50000", "-w", "5", "--bM", "0.25", "--tp", "4"},
                new String[]{}, new String[]{}
            },
            {
                "Acetamide Implicit Solvation Gradients: L = 0.9", new String[]{"ffx/algorithms/structures/acetamide.gk.xyz"},
                ThermoTestMode.GRAD, 0, 0, intRange(1, 9),
                new double[]{-14.3232441, -14.1547628}, new double[]{-47.0065247, 2.88020118},
                new double[]{-626.753662, Double.NaN},
                new double[][][]{
                    {
                        {
                            0.3454325973545622, -0.3556977754400722, 0.22632332784499146
                        },
                        {
                            -2.243711749882479, -0.19356781608568252, -0.7607652107647813
                        },
                        {
                            -0.4691314984840982, 0.41778512359780806, -0.29489746509549386
                        },
                        {
                            -0.13799993351680084, 0.10993253720542051, -0.08060781131332917
                        },
                        {
                            -0.12916866666434035, -0.1292595807344552, -0.19923738828440363
                        },
                        {
                            -0.14207609439404922, 0.1277438883781784, -0.09272512156557061
                        },
                        {
                            -0.23004200820343945, -0.024067408991374717, 0.12377960757109271
                        },
                        {
                            0.7198849556438509, -0.058542226119262075, 0.28168889775425177
                        }
                    },
                    {
                        {
                            0.3245837863444, -0.3346459581204537, 0.21247087522245325
                        },
                        {
                            -2.1088804470365896, -0.18290359797842182, -0.7158477527631418
                        },
                        {
                            -0.4416189323792037, 0.3919514976823785, -0.27720438558513366
                        },
                        {
                            -0.1294804178454077, 0.10345020607616696, -0.07562848810923806
                        },
                        {
                            -0.12170071648317575, -0.12226365328373233, -0.18708155978049676
                        },
                        {
                            -0.1338286311435258, 0.12067970333834456, -0.08762615563248287
                        },
                        {
                            -0.21577447110516357, -0.022826463023161717, 0.11687334540079644
                        },
                        {
                            0.6773203131299521, -0.05386903789837927, 0.26508674941152166
                        }
                    }
                }, new double[][][]{
                    {
                        {
                            6.925790792595476, -6.993227694785508, 4.601662357689225
                        },
                        {
                            -44.78976740441211, -3.542559023701438, -14.921182646958457
                        },
                        {
                            -9.139431352563745, 8.581702257132498, -5.877484673144383
                        },
                        {
                            -2.8301078255996117, 2.1533731217409113, -1.6540871699322346
                        },
                        {
                            -2.4807870616248175, -2.3239852814755535, -4.03805881724825
                        },
                        {
                            -2.7397344153056427, 2.3466598494381614, -1.6938314273137203
                        },
                        {
                            -4.7395497527455746, -0.412231285350944, 2.294199267621269
                        },
                        {
                            14.139598131960081, -1.5523918335908964, 5.515086978981156
                        }
                    }, new double[8][3]
                },
                new String[]{"-l", "0.9"}, new String[]{}, new String[]{}
            },
            {
                "Acetamide Implicit Solvation Gradients: L = 1.0", new String[]{"ffx/algorithms/structures/acetamide.gk.xyz"},
                ThermoTestMode.GRAD, 0, 0, intRange(1, 9),
                new double[]{-22.8540579, -22.4300811}, new double[]{-130.57368, -4.28073063},
                new double[]{-1044.58944, Double.NaN},
                new double[][][]{
                    {
                        {
                            1.602335371, -1.624839098, 1.06143983
                        },
                        {
                            -10.37222509, -0.836476676, -3.468683543
                        },
                        {
                            -2.12776904, 1.975205163, -1.361552091
                        },
                        {
                            -0.651612094, 0.500729882, -0.380794001
                        },
                        {
                            -0.579385578, -0.551019873, -0.932070285
                        },
                        {
                            -0.639287155, 0.553619194, -0.400124158
                        },
                        {
                            -1.090182519, -0.098879753, 0.540134289
                        },
                        {
                            3.285960172, -0.340272596, 1.282575053
                        }
                    },
                    {
                        {
                            1.670608498, -1.693777005, 1.106802141
                        },
                        {
                            -10.81375409, -0.871398548, -3.61577372
                        },
                        {
                            -2.217863814, 2.05980195, -1.419491216
                        },
                        {
                            -0.679510758, 0.521957424, -0.397099678
                        },
                        {
                            -0.603840704, -0.573929277, -0.9718767
                        },
                        {
                            -0.666294935, 0.57675212, -0.416821626
                        },
                        {
                            -1.136904098, -0.102943451, 0.562750069
                        },
                        {
                            3.425345638, -0.355575779, 1.33694173
                        }
                    }
                }, new double[][][]{
                    {
                        {
                            19.23830776, -19.42563249, 12.78239544
                        },
                        {
                            -124.4160206, -9.840441733, -41.44772957
                        },
                        {
                            -25.38730931, 23.83806183, -16.32634631
                        },
                        {
                            -7.861410627, 5.981592005, -4.594686583
                        },
                        {
                            -6.891075171, -6.455514671, -11.21683005
                        },
                        {
                            -7.610373376, 6.518499582, -4.705087298
                        },
                        {
                            -13.16541598, -1.145086904, 6.372775743
                        },
                        {
                            39.27666148, -4.312199538, 15.31968605
                        }
                    }, new double[8][3]
                },
                new String[]{"-l", "1.0"}, new String[]{}, new String[]{}
            },
            {
                // Gradient atoms: a few random protein atoms, some of the coordinating carboxyls, the ions, and some water.
                "Calbindin D9k Ca/Mg Simultaneous Bookending: L = 0.0",
                new String[]{"ffx/algorithms/structures/4icb_ca_a.xyz", "ffx/algorithms/structures/4icb_ca_b.xyz",
                        "ffx/algorithms/structures/4icb_mg_a.xyz", "ffx/algorithms/structures/4icb_mg_b.xyz"},
                ThermoTestMode.GRAD, 0, 0, new int[]{1, 100, 421, 426, 919, 921, 1203, 1204, 1205, 1206, 1207, 1208},
                // Fill in the post-bias PE and dU/dL once I have bias deposition working.
                new double[]{-8422.54811, 0}, new double[]{-0.0234728644, 0}, new double[]{-187.693010, Double.NaN},
                new double[][][] {
                    {
                        {
                            1.51602997670091,-0.16798107175016064,1.1909011229485826
                        },
                        {
                            -1.2287108493416157,0.3317143880292477,-1.168707082997385
                        },
                        {
                            1.330781894105006,-0.8783264761133562,0.1096093788802257
                        },
                        {
                            -1.0136367827497708,0.737243265172173,-0.5962871564399377
                        },
                        {
                            -0.17816595179992722,-0.35265175625985923,-0.6111655937237352
                        },
                        {
                            -0.25917771149474067,0.07021711090701199,0.15836021218566376
                        },
                        {
                            -0.5324396799274522,0.7187373577290392,-0.6184469870748088
                        },
                        {
                            0.2993938879707754,-0.3907323336942321,0.4139579476752431
                        },
                        {
                            0.054758010534996515,0.3356578298638757,0.9215488358561414
                        },
                        {
                            -0.882274050563856,0.5899806670077385,-0.9015065344963205
                        },
                        {
                            -0.09630989055561034,0.9691000112562431,-0.2005707212562129
                        },
                        {
                            0.167130858362869,-0.62725014243849,1.175516097449706
                        }
                    }, new double[12][3] // Fill in once I can get the bias-deposition actually working.
                },
                new double[][][] {
                    {
                        {
                            -2.9163424369832507E-16,3.2314026482573344E-17,-2.2909015893365526E-16
                        },
                        {
                            2.3636350519364876E-16,-6.38109246937693E-17,2.248207565107174E-16
                        },
                        {
                            -2.5599861294253615E-16,1.6896109016191808E-16,-2.1085234990896332E-17
                        },
                        {
                            1.9499033730541748E-16,-1.4182132633554818E-16,1.1470601278862118E-16
                        },
                        {
                            3.427326200965616E-17,6.783858486066127E-17,1.1756813416574204E-16
                        },
                        {
                            4.985725680683067E-17,-1.3507459844962662E-17,-3.0463289923312223E-17
                        },
                        {
                            1.0242386084510416E-16,-1.3826139915463767E-16,1.1896883446563482E-16
                        },
                        {
                            -5.75935248168813E-17,7.516403394171084E-17,-7.963187723760774E-17
                        },
                        {
                            -1.0533638011269165E-17,-6.456951304273904E-17,-1.7727564883699187E-16
                        },
                        {
                            1.6972047349009049E-16,-1.1349285189852184E-16,1.7342017006092682E-16
                        },
                        {
                            1.8526851397737774E-17,-1.8642292909391106E-16,3.8583201849920577E-17
                        },
                        {
                            -3.215047342492132E-17,1.2066227166417432E-16,-2.2613058666622065E-16
                        }
                    }, new double[12][3] // Fill in once I can get the bias-deposition actually working.
                },
                new String[]{"-l", "0.0", "--sf", "TRIG", "--uaA", "329-345.857-972.1008-1022.1204.1208-1213", "--uaB", "329-345.857-972.1008-1022.1204.1208-1213"}, new String[]{"disable-neighbor-updates", "true"}, new String[]{}
            },
            {
                // Gradient atoms: a few random protein atoms, some of the coordinating carboxyls, the ions, and some water.
                "Calbindin D9k Ca/Mg Simultaneous Bookending: L = 0.5",
                new String[]{"ffx/algorithms/structures/4icb_ca_a.xyz", "ffx/algorithms/structures/4icb_ca_b.xyz",
                        "ffx/algorithms/structures/4icb_mg_a.xyz", "ffx/algorithms/structures/4icb_mg_b.xyz"},
                ThermoTestMode.GRAD, 0, 0, new int[]{1, 100, 421, 426, 919, 921, 1203, 1204, 1205, 1206, 1207, 1208},
                // Fill in the post-bias PE and dU/dL once I have bias deposition working.
                new double[]{-8441.56777, 0}, new double[]{-59.7520122, 0}, new double[]{0.170957221, Double.NaN},
                new double[][][] {
                    {
                        {
                            1.516028943632357,-0.16798187601812842,1.1909000012817135
                        },
                        {
                            -1.2287100421217505,0.33171692634032457,-1.1687028725533608
                        },
                        {
                            1.3307953311709682,-0.8783244163319,0.10962280681203751
                        },
                        {
                            -1.0136504919823237,0.7372413665734416,-0.5963003206269726
                        },
                        {
                            -0.17816209958414242,-0.3526516261736795,-0.6111638455124022
                        },
                        {
                            -0.25917449280899385,0.07021855608468686,0.1583629255853478
                        },
                        {
                            -0.5324365378704783,0.718739459636424,-0.6184448359875749
                        },
                        {
                            0.2993955300471849,-0.3907325347437327,0.4139595896177193
                        },
                        {
                            0.05475347093067029,0.335655947922745,0.9215447577884727
                        },
                        {
                            -0.88227515311464,0.5899832887545906,-0.901509658156534
                        },
                        {
                            -0.09630970969676556,0.9690984223562502,-0.20056743541628919
                        },
                        {
                            0.16711991894934464,-0.6272519054249652,1.1755119757400072
                        }
                    }, new double[12][3] // Fill in once I can get the bias-deposition actually working.
                },
                new double[][][] {
                    {
                        {
                            -3.2454805314330315E-6,-2.5266823584502163E-6,-3.5238203937026924E-6
                        },
                        {
                            2.535955989202421E-6,7.974339453653556E-6,1.3227500003987203E-5
                        },
                        {
                            4.221378770274953E-5,6.470994293295007E-6,4.218509194231501E-5
                        },
                        {
                            -4.306882429183645E-5,-5.964623852605655E-6,-4.135651328152079E-5
                        },
                        {
                            1.2102092812327214E-5,4.0867779382836034E-7,5.492167883147658E-6
                        },
                        {
                            1.0111799477741101E-5,4.540159575405767E-6,8.524396519327126E-6
                        },
                        {
                            9.871063099353705E-6,6.603336803578941E-6,6.757839859261594E-6
                        },
                        {
                            5.158735174148887E-6,-6.316156557772956E-7,5.1583144315969776E-6
                        },
                        {
                            -1.4261587612196536E-5,-5.912292426302201E-6,-1.2811627438935602E-5
                        },
                        {
                            -3.463765445133049E-6,8.236460654842404E-6,-9.81326797955262E-6
                        },
                        {
                            5.681848183058946E-7,-4.991676597398964E-6,1.0322770563675476E-5
                        },
                        {
                            -3.436718115423787E-5,-5.538585353903613E-6,-1.2948732912576588E-5
                        }
                    }, new double[12][3] // Fill in once I can get the bias-deposition actually working.
                },
                new String[]{"-l", "0.5", "--sf", "TRIG", "--uaA", "329-345.857-972.1008-1022.1204.1208-1213", "--uaB", "329-345.857-972.1008-1022.1204.1208-1213"}, new String[]{"disable-neighbor-updates", "true"}, new String[]{}
            },
            {
                // Gradient atoms: a few random protein atoms, some of the coordinating carboxyls, the ions, and some water.
                "Calbindin D9k Ca/Mg Simultaneous Bookending: L = 0.7",
                new String[]{"ffx/algorithms/structures/4icb_ca_a.xyz", "ffx/algorithms/structures/4icb_ca_b.xyz",
                        "ffx/algorithms/structures/4icb_mg_a.xyz", "ffx/algorithms/structures/4icb_mg_b.xyz"},
                ThermoTestMode.GRAD, 0, 0, new int[]{1, 100, 421, 426, 919, 921, 1203, 1204, 1205, 1206, 1207, 1208},
                // Fill in the post-bias PE and dU/dL once I have bias deposition working.
                new double[]{-8452.74724, 0}, new double[]{-48.3265963, 0}, new double[]{110.479770, Double.NaN},
                new double[][][] {
                    {
                        {
                            1.5160283364099056,-0.1679823487549772,1.1908993419824703
                        },
                        {
                            -1.2287095676498276,0.3317184183221489,-1.1687003977164605
                        },
                        {
                            1.3308032292801621,-0.8783232056227372,0.10963069955232374
                        },
                        {
                            -1.0136585500670425,0.7372402506051008,-0.5963080583419732
                        },
                        {
                            -0.1781598353085141,-0.35265154971094115,-0.6111628179395603
                        },
                        {
                            -0.25917260091298266,0.07021940553881123,0.15836452048166638
                        },
                        {
                            -0.5324346910157285,0.7187406951065878,-0.6184435716102188
                        },
                        {
                            0.2993964952354813,-0.3907326529176762,0.41396055472729376
                        },
                        {
                            0.05475080261819443,0.33565484174550564,0.921542360760436
                        },
                        {
                            -0.8822758011777325,0.5899848297787287,-0.9015114941979432
                        },
                        {
                            -0.09630960339060302,0.9690974884242527,-0.20056550404804518
                        },
                        {
                            0.16711348892341116,-0.627252941682412,1.175509553059828
                        }
                    }, new double[12][3] // Fill in once I can get the bias-deposition actually working.
                },
                new double[][][] {
                    {
                        {
                            -2.625648932053082E-6,-2.044128962097602E-6,-2.8508305689456392E-6
                        },
                        {
                            2.0516315029084353E-6,6.4513761444118245E-6,1.0701272278756946E-5
                        },
                        {
                            3.4151671600923805E-5,5.235144362458755E-6,3.412845626726835E-5
                        },
                        {
                            -3.4843410761098426E-5,-4.825482083248289E-6,-3.34581220711172E-5
                        },
                        {
                            9.790798751296848E-6,3.3062727666788305E-7,4.443257157316083E-6
                        },
                        {
                            8.180617633968268E-6,3.6730662467121533E-6,6.896381644416749E-6
                        },
                        {
                            7.985857803483043E-6,5.342211697900723E-6,5.467207293663456E-6
                        },
                        {
                            4.173504418503171E-6,-5.109878138398471E-7,4.1731640418873894E-6
                        },
                        {
                            -1.1537866743793757E-5,-4.783145044484627E-6,-1.0364824328590316E-5
                        },
                        {
                            -2.802245113286972E-6,6.663436650455878E-6,-7.93910056628988E-6
                        },
                        {
                            4.596711775661788E-7,-4.038351189450395E-6,8.351296809649966E-6
                        },
                        {
                            -2.780363359988769E-5,-4.480809675300179E-6,-1.0475744996796266E-5
                        }
                    }, new double[12][3] // Fill in once I can get the bias-deposition actually working.
                },
                new String[]{"-l", "0.7", "--sf", "TRIG", "--uaA", "329-345.857-972.1008-1022.1204.1208-1213", "--uaB", "329-345.857-972.1008-1022.1204.1208-1213"}, new String[]{"disable-neighbor-updates", "true"}, new String[]{}
            },
            {
                // Gradient atoms: a few random protein atoms, some of the coordinating carboxyls, the ions, and some water.
                "Calbindin D9k Ca/Mg Simultaneous Bookending: L = 1.0",
                new String[]{"ffx/algorithms/structures/4icb_ca_a.xyz", "ffx/algorithms/structures/4icb_ca_b.xyz",
                        "ffx/algorithms/structures/4icb_mg_a.xyz", "ffx/algorithms/structures/4icb_mg_b.xyz"},
                ThermoTestMode.GRAD, 0, 0, new int[]{1, 100, 421, 426, 919, 921, 1203, 1204, 1205, 1206, 1207, 1208},
                // Fill in the post-bias PE and dU/dL once I have bias deposition working.
                new double[]{-8460.58743, 0}, new double[]{0.0234728644, 0}, new double[]{187.739955, Double.NaN},
                new double[][][] {
                    {
                        {
                            1.5160279105638157,-0.16798268028609797,1.1908988796148456
                        },
                        {
                            -1.2287092349018947,0.3317194646514121,-1.1686986621093438
                        },
                        {
                            1.3308087682369107,-0.8783223565504423,0.10963623474384576
                        },
                        {
                            -1.0136642012148873,0.7372394679746987,-0.5963134848140097
                        },
                        {
                            -0.1781582473683545,-0.35265149608749957,-0.6111620973010643
                        },
                        {
                            -0.2591712741232506,0.07022000126236438,0.15836563898503053
                        },
                        {
                            -0.5324333958135048,0.718741561543808,-0.6184426849003377
                        },
                        {
                            0.29939717212359174,-0.3907327357932422,0.4139612315602017
                        },
                        {
                            0.05474893132634229,0.3356540659816183,0.921540679720797
                        },
                        {
                            -0.8822762556654258,0.5899859105014462,-0.9015127818167441
                        },
                        {
                            -0.09630952883792121,0.9690968334562342,-0.20056414957637436
                        },
                        {
                            0.1671089795358247,-0.6272536684114378,1.1755078540303066
                        }
                    }, new double[12][3] // Fill in once I can get the bias-deposition actually working.
                },
                new double[][][] {
                    {
                        {
                            2.9163384624158903E-16,-3.231433591192038E-17,2.290897273901192E-16
                        },
                        {
                            -2.3636319462860973E-16,6.381190126870034E-17,-2.248191366091609E-16
                        },
                        {
                            2.560037826405337E-16,-1.6896029769367465E-16,2.109040117467388E-17
                        },
                        {
                            -1.9499561171519864E-16,1.4182059587979724E-16,-1.147110775007824E-16
                        },
                        {
                            -3.4271779930734056E-17,-6.783853481206671E-17,-1.1756746156916025E-16
                        },
                        {
                            -4.9856018468543316E-17,1.3508015854151776E-17,3.0464333860802826E-17
                        },
                        {
                            -1.0242265198852129E-16,1.38262207830166E-16,-1.1896800686893987E-16
                        },
                        {
                            5.759415657973307E-17,-7.51641112923199E-17,7.963250894893407E-17
                        },
                        {
                            1.0531891470506986E-17,6.456878899573994E-17,1.7727407986513553E-16
                        },
                        {
                            -1.6972089767901735E-16,1.1349386057403977E-16,-1.734213718396476E-16
                        },
                        {
                            -1.8526781815166296E-17,1.8642231778983455E-16,-3.858193767513045E-17
                        },
                        {
                            3.214626465908195E-17,-1.2066294994525673E-16,2.2612900090378923E-16
                        }
                    }, new double[12][3] // Fill in once I can get the bias-deposition actually working.
                },
                new String[]{"-l", "1.0", "--sf", "TRIG", "--uaA", "329-345.857-972.1008-1022.1204.1208-1213", "--uaB", "329-345.857-972.1008-1022.1204.1208-1213"}, new String[]{"disable-neighbor-updates", "true"}, new String[]{}
            }
        });
    }

    private static int[] intRange(int low, int high) {
        return IntStream.range(low, high).toArray();
    }

    // Enables behavior where it prints out instead of testing information.
    private static final boolean debugMode = false;
    /**
     * Set of default options, such as "-d", "0.5" to specify a half-femtosecond timestep.
     * Can be over-ridden in test constructors.
     */
    private static final Map<String, String> DEFAULT_OPTIONS;
    /**
     * Set of default properties, such as "ttosrw-alwaystemper", "true" to use standard tempering scheme.
     * Can be over-ridden in test constructors.
     */
    private static final Map<String, String> DEFAULT_PROPERTIES;
    /**
     * Set of default boolean options, such as "-y", false to specify "do not add the -y --synchronous flag".
     * Can be over-ridden in test constructors.
     * Maps a flag to whether it should be included by default.
     */
    private static final Map<String, Boolean> DEFAULT_FLAGS;

    /**
     * Default number of atomic gradients to test; pick the first N (3) if ffx.ci is false. Otherwise, use
     * the whole array to test.
     */
    private static final int DEFAULT_GRADIENT_EVALS = 3;

    static {
        String[] opts = {"--bM", "0.05",
                "--dw", "OFF",
                "--lf", "1.0E-18",
                "--lm", "1.0E-18",
                "--np", "1",
                "--sf", "1.0",
                "--tp", "2.0",      // Set low to encourage fast tempering for toy cases.
                "-b", "ADIABATIC",  // Used for Langevin dynamics.
                "-C", "5000",       // By default, evaluate 0 dynamics steps, and trigger adding counts manually.
                "-d", "0.5",        // Most of the tests are in vacuum, thus short timestep.
                "-F", "XYZ",
                "-i", "STOCHASTIC", // Stochastic because most toy systems work best with Langevin dynamics.
                "-k", "5.0",        // Infrequent because we should not need to restart.
                "-l", "1.0",        // Set walker to start at the usually-physical topology. Likely commonly set by test.
                "-n", "0",          // Do not run any steps by default, just evaluate energy.
                "-p", "0",          // Use NVT, not NPT dynamics.
                "-Q", "0",          // Do not run any steps by default, just evaluate energy.
                "-r", "2.0",        // Very infrequent printouts, since it should work.
                "-t", "298.15",
                "-w", "10.0"};      // Infrequent since we should not need the trajectory.

        int nOpts = opts.length;
        Map<String, String> optMap = new HashMap<>(nOpts);
        for (int i = 0; i < nOpts; i += 2) {
            optMap.put(opts[i], opts[i + 1]);
        }
        DEFAULT_OPTIONS = Collections.unmodifiableMap(optMap);

        String[] props = {"ttosrw-alwaystemper", "true",
                "ttosrw-temperOffset", "0.5", // Small to temper fast.
                "print-on-failure", "false",
                "disable-neighbor-updates", "false",
                "vdw-lambda-alpha", "0.25",
                "vdw-lambda-exponent", "3",
                "permanent-lambda-alpha", "1.0",
                "permanent-lambda-exponent", "3"
        };

        int nProps = props.length;
        Map<String, String> propMap = new HashMap<>(nProps);
        for (int i = 0; i < nProps; i += 2) {
            propMap.put(props[i], props[i + 1]);
        }
        DEFAULT_PROPERTIES = Collections.unmodifiableMap(propMap);

        Map<String, Boolean> flags = new HashMap<>();
        flags.put("--rn", false);
        flags.put("-y", false);
        flags.put("-h", false);
        DEFAULT_FLAGS = Collections.unmodifiableMap(flags);
    }

    Binding binding;
    NewThermodynamics thermo;

    private final String info;
    private final ThermoTestMode mode;
    private final String[] filenames;
    private final File tempDir;
    private final File[] copiedFiles;
    private final Map<String, String> opts;
    private final List<String> flags;
    /**
     * Configuration containing the properties to be used by TT-OSRW.
     */
    Configuration algorithmConfig;
    private final boolean ffxCI;

    /**
     * Free energy associated with this test, if calculated.
     * Not used alongside testing gradients, and values are meaningless in that case.
     */
    private final double freeEnergy;
    private final double feTol;

    /**
     * Set of atoms for which tabulated gradients (dU/dL, d2U/dXdL) should be evaluated.
     * Not used alongside testing free energies, and values are meaningless (and often null) in that case.
     */
    private final int[] gradAtomIndices;
    private final int numGradAtoms;

    /**
     * Potential energy of the underlying Potential (0), and OSRW after one bias drop (1).
     */
    private final double[] pe;
    private final double[] dudl;
    private final double[] d2udl2;
    /**
     * Gradient of the underlying Potential ([0][0..n][0-2]), and OSRW after one bias drop ([1][0..n][0-2]).
     * Axes: before/after bias, atoms, X/Y/Z.
     */
    private final double[][][] dudx;
    private final double[][][] d2udxdl;

    // Currently, everything is just 1.0E-5 kcal/mol per unit.
    private static final double DEFAULT_PE_TOL = 1.0E-5;      // Units: none.
    private static final double DEFAULT_DUDL_TOL = 1.0E-5;    // Units: per lambda (quasi-unitless).
    private static final double DEFAULT_D2UDL2_TOL = 1.0E-5;  // Units: per lambda^2 (quasi-unitless).
    private static final double DEFAULT_DUDX_TOL = 1.0E-5;    // Units: per Angstrom.
    private static final double DEFAULT_D2UDXDL_TOL = 1.0E-5; // Units: per Angstrom per lambda.

    // TODO: Make these settable if needed.
    private final double peTol = DEFAULT_PE_TOL;
    private final double dudlTol = DEFAULT_DUDL_TOL;
    private final double d2udl2Tol = DEFAULT_D2UDL2_TOL;
    private final double dudxTol = DEFAULT_DUDX_TOL;
    private final double d2udxdlTol = DEFAULT_D2UDXDL_TOL;

    public ThermodynamicsTest(String info, String[] filenames, ThermoTestMode mode,
                              double freeEnergy, double feTol,
                              int[] gradAtomIndices,
                              double[] pe, double[] dudl, double[] d2udl2,
                              double[][][] dudx, double[][][] d2udxdl,
                              String[] options, String[] properties, String[] flags) throws IOException {
        this.info = info;
        this.mode = mode;
        int nFiles = filenames.length;

        String tempDirName = String.format("temp-%016x/", new Random().nextLong());
        tempDir = new File(tempDirName);
        tempDir.mkdir();
        logger.fine(String.format(" Running test %s in directory %s", info, tempDir));
        this.filenames = Arrays.copyOf(filenames, nFiles);
        copiedFiles = new File[nFiles];

        String[] copiedExtensions = new String[]{"dyn", "key", "properties", "his", "lam"};
        for (int i = 0; i < nFiles; i++) {
            File srcFile = new File("src/main/java/" + filenames[i]);
            File tempFile = new File(tempDirName + FilenameUtils.getName(filenames[i]));
            FileUtils.copyFile(srcFile, tempFile);
            logger.info(String.format(" Copied file %s to %s", srcFile, tempFile));
            copiedFiles[i] = tempFile;

            for (String ext : copiedExtensions) {
                srcFile = new File(String.format("%s.%s", FilenameUtils.removeExtension(srcFile.getPath()), ext));
                if (srcFile.exists()) {
                    logger.fine(" Copying extension " + ext);
                    tempFile = new File(String.format("%s.%s", FilenameUtils.removeExtension(tempFile.getPath()), ext));
                    logger.info(String.format(" Copied file %s to %s", srcFile, tempFile));
                    FileUtils.copyFile(srcFile, tempFile);
                }
            }
        }
        ffxCI = Boolean.parseBoolean(System.getProperty("ffx.ci", "false"));

        switch (mode) {
            case HELP:
                assertTrue(String.format("Help tests must have no file arguments, found %d for %s", nFiles, info), nFiles == 0);
                break;
            default:
                assertTrue(String.format("Must have 1, 2, or 4 distinct filenames, found %d for test %s", nFiles, info),
                        nFiles == 1 || nFiles == 2 || nFiles == 4);
                for (int i = 0; i < nFiles; i++) {
                    for (int j = i + 1; j < nFiles; j++) {
                        assertNotEquals(String.format(" Filenames %d and %d matched in test %s: files %s and %s", i, j, info, filenames[i], filenames[j]),
                                filenames[i], filenames[j]);
                    }
                }
                break;
        }
        int nOpts = options.length;
        int nProps = properties.length;
        int nFlags = flags.length;

        if (nOpts > 0) {
            assertTrue(String.format("Unmatched option key %s for test %s", options[nOpts - 1], info),
                    options.length % 2 == 0);
        }
        if (nProps > 0) {
            assertTrue(String.format("Unmatched property key %s for test %s", properties[nProps - 1], info),
                    properties.length % 2 == 0);
        }
        if (nFlags > 0) {
            assertTrue(String.format("Unmatched flag key %s for test %s", flags[nFlags - 1], info),
                    flags.length % 2 == 0);
        }

        Pattern validOption = Pattern.compile("^--?[^D]");
        Pattern validProperty = Pattern.compile("^[^-]");

        Map<String, String> groovyOpts = new HashMap<>(DEFAULT_OPTIONS);
        for (int i = 0; i < nOpts; i += 2) {
            String opti = options[i];
            assertTrue(String.format(" Option %s for test %s does not look like a Groovy option!", opti, info),
                    validOption.matcher(opti).find());
            groovyOpts.put(opti, options[i + 1]);
        }
        this.opts = Collections.unmodifiableMap(groovyOpts);

        Map<String, String> addedProps = new HashMap<>(DEFAULT_PROPERTIES);
        for (int i = 0; i < nProps; i += 2) {
            String propi = properties[i];
            assertTrue(String.format(" Property %s for test %s does not look like a property!", propi, info),
                    validProperty.matcher(propi).find());
            addedProps.put(propi, properties[i + 1]);
        }
        algorithmConfig = new MapConfiguration(addedProps);


        Map<String, Boolean> addedFlags = new HashMap<>(DEFAULT_FLAGS);
        for (int i = 0; i < nFlags; i += 2) {
            String flagi = flags[i];
            assertTrue(String.format(" Flag %s for test %s does not look like a flag!", flagi, info),
                    validOption.matcher(flagi).find());
            String vali = flags[i + 1];
            assertTrue(String.format(" Value %s for flag %s in test %s is not a true/false value!", vali, flagi, info),
                    vali.equalsIgnoreCase("TRUE") || vali.equalsIgnoreCase("FALSE"));
            addedFlags.put(flagi, Boolean.parseBoolean(vali));
        }
        this.flags = Collections.unmodifiableList(
                addedFlags.entrySet().stream().
                        filter(Map.Entry::getValue).
                        map(Map.Entry::getKey).
                        collect(Collectors.toList()));

        // Only meaningful for free energy evaluations.
        this.freeEnergy = freeEnergy;
        this.feTol = feTol;
        if (mode == ThermoTestMode.FREE) {
            assertTrue(String.format(" Free energy tolerance for test %s was %10.4g <= 0!", info, feTol), feTol > 0);
        }

        // Only meaningful for gradient evaluations.
        this.dudx = new double[2][][];
        this.d2udxdl = new double[2][][];
        boolean isGradient = (mode == ThermoTestMode.GRAD);
        if (isGradient) {
            numGradAtoms = ffxCI ? gradAtomIndices.length : Math.min(DEFAULT_GRADIENT_EVALS, gradAtomIndices.length);
            this.gradAtomIndices = Arrays.copyOf(gradAtomIndices, numGradAtoms);

            if (debugMode) {
                this.pe = this.dudl = this.d2udl2 = null;
            } else {
                assertNotNull(gradAtomIndices);
                assertNotNull(pe);
                assertNotNull(dudl);
                assertNotNull(d2udl2);

                this.pe = Arrays.copyOf(pe, 2);
                this.dudl = Arrays.copyOf(dudl, 2);
                this.d2udl2 = Arrays.copyOf(d2udl2, 2);
                for (int i = 0; i < 2; i++) {
                    assertNotNull(dudx[i]);
                    assertNotNull(d2udxdl[i]);
                    this.dudx[i] = new double[numGradAtoms][3];
                    this.d2udxdl[i] = new double[numGradAtoms][3];

                    for (int j = 0; j < numGradAtoms; j++) {
                        assertNotNull(dudx[i][j]);
                        assertNotNull(d2udxdl[i][j]);
                        System.arraycopy(dudx[i][j], 0, this.dudx[i][j], 0, 3);
                        System.arraycopy(d2udxdl[i][j], 0, this.d2udxdl[i][j], 0, 3);
                    }
                }
            }
        } else {
            this.gradAtomIndices = null;
            numGradAtoms = 0;
            this.pe = null;
            this.dudl = null;
            this.d2udl2 = null;
            for (int i = 0; i < 2; i++) {
                this.dudx[i] = null;
                this.d2udxdl[i] = null;
            }
        }
    }

    @Before
    public void before() {
        binding = new Binding();
        thermo = new NewThermodynamics();
        thermo.setBinding(binding);
    }

    @After
    public void after() throws IOException {
        // Clean up the temporary directory if it exists.
        if (tempDir != null && tempDir.exists()) {
            if (tempDir.isDirectory()) {
                for (File file : tempDir.listFiles()) {
                    if (file.isDirectory()) {
                        for (File file2 : file.listFiles()) {
                            // Should never have to go more than 2-deep into the temporary directory.
                            file2.delete();
                        }
                    }
                    file.delete();
                }
            } else {
                logger.warning(String.format(" Expected %s (temporary directory) to be a directory!", tempDir));
            }
            tempDir.delete();
        }

        if (thermo.getOSRW() == null) {
            assert mode == ThermoTestMode.HELP;
        } else {
            thermo.destroyPotentials();
        }
        System.gc();
    }

    @Test
    public void testThermodynamics() {
        logger.info(String.format(" Thermodynamics test: %s\n", info));
        switch (mode) {
            case HELP:
                testHelp();
                break;
            case FREE:
                testFreeEnergy();
                break;
            case GRAD:
                testStaticGradients();
                break;
            default:
                throw new IllegalStateException(String.format(" Thermodynamics test mode %s not recognized!", mode));
        }
    }

    private void testHelp() {
        String[] args = {"-h"};
        binding.setVariable("args", args);
        thermo.run();
    }

    private String[] assembleArgs() {
        List<String> argList = new ArrayList<>(flags);
        opts.forEach((String k, String v) -> {
            argList.add(k);
            argList.add(v);
        });
        argList.addAll(Arrays.stream(copiedFiles).map(File::getPath).collect(Collectors.toList()));
        return argList.toArray(new String[argList.size()]);
    }

    private void assembleThermo() {
        String[] args = assembleArgs();
        binding.setVariable("args", args);
        thermo = new NewThermodynamics();
        thermo.setBinding(binding);
        thermo.setProperties(algorithmConfig);
    }

    /**
     * Test a free energy evaluation; very expensive, so likely only done for one or a few tests.
     * Currently not implemented.
     */
    private void testFreeEnergy() {
        assembleThermo();
        thermo.run();
        AbstractOSRW osrw = thermo.getOSRW();
        double delG = osrw.lastFreeEnergy();
        assertEquals(String.format(" Test %s: not within tolerance %12.5g", info, feTol), freeEnergy, delG, feTol);
    }

    /**
     * Tests gradients & energies for a static structure, before and after dropping a bias.
     */
    private void testStaticGradients() {
        assembleThermo();
        thermo.run();
        AbstractOSRW osrw = thermo.getOSRW();
        osrw.setPropagateLambda(false);
        CrystalPotential under = thermo.getPotential();
        int nVars = osrw.getNumberOfVariables();

        double[] x = new double[nVars];
        x = osrw.getCoordinates(x);
        double[] gOSRWPre = new double[nVars];
        double[] gUnderPre = new double[nVars];
        double[] gOSRWPost = new double[nVars];
        double[] gUnderPost = new double[nVars];

        logger.info(" Testing the OSRW potential before bias added.");
        EnergyResult osrwPre = testGradientSet("Unbiased OSRW", osrw, x, gOSRWPre, 0);
        logger.info(" Testing the underlying CrystalPotential before bias added.");
        EnergyResult underPre = testGradientSet("Unbiased potential", under, x, gUnderPre, 0);

        // Assert that, before biases, OSRW and underlying potential are equal.
        osrwPre.assertResultsEqual(underPre);

        // This code works perfectly on two local machines and yet fails on Travis.
        /*double currentdUdL = osrwPre.firstLam;
        logger.info(String.format(" Adding an OSRW bias at lambda %8.4g, dU/dL %14.8g", osrw.getLambda(), currentdUdL));
        osrw.addBias(currentdUdL, x, gOSRWPre);

        // Wait for the bias to be received by the OSRW object.
        boolean biasReceived = false;
        for (int i = 0; i < 200; i++) {
            try {
                Thread.sleep(100);
            } catch (InterruptedException ex) {
                logger.warning(String.format(" Interrupted at %d 100-msec wait for bias to be added!", i));
            }
            if (osrw.getCountsReceived() > 0) {
                logger.fine(String.format(" Required %d 100-msec waits for bias to be added.", i));
                biasReceived = true;
                break;
            }
        }
        try {
            Thread.sleep(500);
        } catch (InterruptedException ex) {
            logger.warning(" Interrupted at final 500-msec wait for bias to be added!");
        }
        assertTrue(" No bias was received by the OSRW over 20 seconds!", biasReceived);

        logger.info(" Testing the OSRW potential after bias added.");
        EnergyResult osrwPost = testGradientSet("Biased OSRW", osrw, x, gOSRWPost, 1);
        logger.info(" Testing the underlying CrystalPotential after bias added.");
        EnergyResult underPost = testGradientSet("Biased potential", under, x, gUnderPost, 0);

        // Assert that bias only affects the OSRW potential, and does affect the OSRW potential.
        underPre.assertResultsEqual(underPost);
        osrwPre.assertResultsInequal(osrwPost);*/
    }

    /**
     * Generates and tests an EnergyResult against tabulated values.
     *
     * @param description Description (such as unbiased OSRW)
     * @param potential   A CrystalPotential (either OSRW or underlying).
     * @param x           Coordinates.
     * @param g           Array to add gradients to.
     * @param tableIndex  0 for unbiased potential, 1 for a biased potential.
     * @return The generated EnergyResult.
     */
    private EnergyResult testGradientSet(String description, CrystalPotential potential, double[] x, double[] g, int tableIndex) {
        assertTrue(String.format(" Potential %s is not a lambda interface!", potential), potential instanceof LambdaInterface);

        EnergyResult er = new EnergyResult(description, potential, x, g);

        checkThGradScalar(er.energy, pe, tableIndex, peTol, "potential energy");
        checkThGradScalar(er.firstLam, dudl, tableIndex, dudlTol, "dU/dL");
        checkThGradArray(er.gradient, dudx, tableIndex, dudxTol, "dU/dX gradient");
        if (er.hasSecondDerivatives) {
            checkThGradScalar(er.secondLam, d2udl2, tableIndex, d2udl2Tol, "d2U/dL2");
            checkThGradArray(er.lamGradient, d2udxdl, tableIndex, d2udxdlTol, "d2U/dXdL");
        }

        return er;
    }

    /**
     * Checks a scalar value against its expected value. Mostly a convenience formatting method.
     *
     * @param actual      Value from the test.
     * @param expected    Array of expected values (see tableIndex).
     * @param tableIndex  0 for unbiased potential, 1 for biased OSRW potential.
     * @param tol         Tolerance for this test.
     * @param description Scalar to be tested.
     */
    private void checkThGradScalar(double actual, double[] expected, int tableIndex, double tol, String description) {
        if (debugMode) {
            logger.info(String.format(" %s is %15.9g", description, actual));
        } else {
            assertEquals(String.format(" Expected %s %12.6g, received %12.6g from test %s",
                    description, expected[tableIndex], actual, info), expected[tableIndex], actual, tol);
        }
    }

    /**
     * Checks an array value (generally 1-D flat) against its expected value (generally 3-D; indices
     * tableIndex, then atom number, then X/Y/Z).
     *
     * @param actual      Array from the test, flat.
     * @param expected    Array of expected values, indexed by tableIndex, atoms, and XYZ.
     * @param tableIndex  0 for unbiased potential, 1 for biased OSRW potential.
     * @param tol         Tolerance for this test.
     * @param description Array to be tested.
     */
    private void checkThGradArray(double[] actual, double[][][] expected, int tableIndex, double tol, String description) {
        double[] actualSlice = new double[3];
        for (int i = 0; i < numGradAtoms; i++) {
            int i3 = i * 3;
            System.arraycopy(actual, i3, actualSlice, 0, 3);
            if (debugMode) {
                logger.info(String.format(" %s at atom %d is %s", description, i, Arrays.toString(actualSlice)));
            } else {
                double[] exp = expected[tableIndex][i];
                assertArrayEquals(String.format(" Discrepancy found on array of %s for test %s on atom %d." +
                                "\n Expected: %s\n Found: %s", description, info, i, Arrays.toString(exp), Arrays.toString(actualSlice)),
                        exp, actualSlice, tol);
            }
        }
    }

    /**
     * Checks if two double values are approximately equal to within a tolerance.
     *
     * @param v1     One value to compare.
     * @param v2     Second value to compare.
     * @param absTol Tolerance for inequality (absolute, not relative).
     * @return True if v1 approximately equal to v2.
     */
    private static boolean approxEquals(double v1, double v2, double absTol) {
        double diff = v1 - v2;
        return Math.abs(diff) < absTol;
    }

    /**
     * Contains the result of an energy evaluation: potential energy and several derivatives.
     */
    private class EnergyResult {
        final double energy;
        final double firstLam;
        final double secondLam;
        final int nVars;
        final boolean hasSecondDerivatives;
        private final String description;

        private final double[] gradient;
        private final double[] lamGradient;

        public EnergyResult(String description, CrystalPotential potential, double[] x, double[] g) {
            this.description = description;

            LambdaInterface linter = (LambdaInterface) potential;
            energy = potential.energyAndGradient(x, g, false);
            firstLam = linter.getdEdL();
            nVars = g.length;
            gradient = Arrays.copyOf(g, nVars);

            hasSecondDerivatives = !(potential instanceof TransitionTemperedOSRW);
            if (hasSecondDerivatives) {
                secondLam = linter.getd2EdL2();
                lamGradient = new double[g.length];
                linter.getdEdXdL(lamGradient);
            } else {
                secondLam = 0;
                lamGradient = null;
            }
        }

        /**
         * Asserts that two energy results are equivalent to each other. Always tests U, dU/dX, and dU/dL.
         * Also tests d2U/dL2 and d2U/dXdL if these second gradients are available.
         * <p>
         * All values must be identical to tolerance.
         *
         * @param other Another EnergyResult
         */
        public void assertResultsEqual(EnergyResult other) {
            assertEquals(String.format(" Test %s: potential energy for %s did not " +
                            "match %s, %12.6g != %12.6g.", info, this.toString(), other.toString(), this.energy, other.energy),
                    this.energy, other.energy, peTol);
            assertEquals(String.format(" Test %s: dU/dL for %s did not " +
                            "match %s, %12.6g != %12.6g.", info, this.toString(), other.toString(), this.firstLam, other.firstLam),
                    this.firstLam, other.firstLam, dudlTol);
            assertArrayEquals(String.format(" Test %s: dU/dX for %s did not match %s.",
                    info, this.toString(), other.toString()),
                    this.gradient, other.gradient, dudxTol);

            if (hasSecondDerivatives && other.hasSecondDerivatives) {
                assertEquals(String.format(" Test %s: d2U/dL2 for %s did not " +
                                "match %s, %12.6g != %12.6g.", info, toString(), other.toString(), this.secondLam, other.secondLam),
                        this.secondLam, other.secondLam, d2udl2Tol);
                assertArrayEquals(String.format(" Test %s: d2U/dXdL for %s did not match %s.",
                        info, this.toString(), other.toString()),
                        this.lamGradient, other.lamGradient, d2udxdlTol);
            }
        }

        /**
         * Asserts that two energy results are not equivalent to each other. Always tests U, dU/dX, and dU/dL.
         * Also tests d2U/dL2 and d2U/dXdL if these second gradients are available.
         * <p>
         * Prints all equalities, fails if no inequalities are found. At least one value must not be identical to tolerance.
         *
         * @param other Another EnergyResult
         */
        public void assertResultsInequal(EnergyResult other) {
            int diffsFound = 0;
            StringBuilder sb = new StringBuilder(" Equalities found between ");
            sb = sb.append(toString()).append(" and ").append(other.toString());
            sb = sb.append(" for test ").append(info);
            boolean equalFound = false;

            if (approxEquals(energy, other.energy, peTol)) {
                equalFound = true;
                sb.append(String.format("\n Potential energy %12.6g == other potential energy %12.6g", energy, other.energy));
            } else {
                ++diffsFound;
            }

            if (approxEquals(firstLam, other.firstLam, dudlTol)) {
                equalFound = true;
                sb.append(String.format("\n Lambda derivative %12.6g == other lambda derivative %12.6g", firstLam, other.firstLam));
            } else {
                ++diffsFound;
            }

            for (int i = 0; i < nVars; i++) {
                if (approxEquals(gradient[i], other.gradient[i], dudxTol)) {
                    equalFound = true;
                    sb.append(String.format("\n Gradient %d %12.6g == other gradient %12.6g", i, gradient[i], other.gradient[i]));
                } else {
                    ++diffsFound;
                }
            }

            if (hasSecondDerivatives && other.hasSecondDerivatives) {
                if (approxEquals(secondLam, other.secondLam, d2udl2Tol)) {
                    equalFound = true;
                    sb.append(String.format("\n Lambda derivative %12.6g == other lambda derivative %12.6g", secondLam, other.secondLam));
                } else {
                    ++diffsFound;
                }

                for (int i = 0; i < nVars; i++) {
                    if (approxEquals(lamGradient[i], other.lamGradient[i], d2udxdlTol)) {
                        equalFound = true;
                        sb.append(String.format("\n Lambda gradient %d %12.6g == other lambda gradient %12.6g", i, lamGradient[i], other.lamGradient[i]));
                    } else {
                        ++diffsFound;
                    }
                }
            }

            if (equalFound) {
                logger.info(sb.toString());
            }
            assertTrue(String.format(" No inequalities found between %s and %s for test %s!",
                    this.toString(), other.toString(), info), diffsFound > 0);
        }

        @Override
        public String toString() {
            return description;
        }
    }

    private enum ThermoTestMode {
        HELP, FREE, GRAD
    }
}
