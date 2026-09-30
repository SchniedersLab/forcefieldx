//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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
package ffx.potential.commands;

import ffx.potential.bonded.Atom;
import ffx.potential.cli.PotentialCommand;
import ffx.utilities.FFXBinding;

import java.util.Locale;
import java.util.Scanner;
import java.util.logging.Level;

import org.beehive.jitllm.Options;
import org.beehive.jitllm.api.FinishReason;
import org.beehive.jitllm.api.GenerationRequest;
import org.beehive.jitllm.api.GenerationResult;
import org.beehive.jitllm.api.GenerationSession;
import org.beehive.jitllm.api.LocalModel;
import org.beehive.jitllm.api.LocalModels;
import org.beehive.jitllm.api.ModelOptions;
import org.beehive.jitllm.api.TextGenerationModel;
import org.beehive.jitllm.auxiliary.RunMetrics;
import org.beehive.jitllm.integration.cli.ModelRunConfig;
import org.beehive.jitllm.integration.cli.StartupDiagnostics;

import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

/**
 * Chat with an LLM about this molecule.
 * <p>
 * Usage:
 * ffxc Chat [options] &lt;message&gt;
 */
@Command(name = "Chat", description = "Chat with an LLM.")
public class Chat extends PotentialCommand {

  public static final boolean SHOW_PERF_INTERACTIVE =
      Boolean.parseBoolean(
          System.getProperty(
              "jitllm.ShowPerfInteractive",
              "false")); // Show performance metrics in interactive mode

  /*
        out.println("Usage: jitllm run|chat [options] (legacy mode flags remain supported)");
        out.println();
        out.println("Options:");
        out.println("  --model, -m <path>            required, path to .gguf file");
        out.println("  --interactive, --chat, -i     run in chat mode");
        out.println("  --instruct                    run in instruct (once) mode, default mode");
        out.println("  --prompt, -p <string>         input prompt");
        out.println("  --system-prompt, -sp <string> (optional) system prompt (Llama models)");
        out.println("  --suffix <string>             suffix for fill-in-the-middle request (Codestral)");
        out.println("  --temperature, -temp <float>  temperature in [0,inf], default: auto-detected from model family");
        out.println("  --top-p <float>               p value in top-p (nucleus) sampling in [0,1], default: auto-detected from model family");
        out.println("  --seed <long>                 random seed, default System.nanoTime()");
        out.println("  --ctx-size, -c <int>          context capacity, including prompt and generated tokens");
        out.println("  --max-new-tokens <int>        generated-token limit per request/turn, default: context capacity");
        out.println("  --max-tokens, -n <int>        capacity in positions, prompt plus generated tokens (the context"
                        + " this run allocates; a prompt this long or longer is refused), default "
                        + DEFAULT_MAX_TOKENS);
        out.println("  --stream <boolean>            print tokens during generation; may cause encoding artifacts for non ASCII text, default true");
        out.println("  --echo <boolean>              print ALL tokens to stderr, if true, recommended to set --stream=false, default false");
        out.println("  --with-prefill-decode         enable prefill/decode separation (skip logits during prefill)");
        out.println("  --batch-prefill-size <int>    batched prefill chunk size; requires --with-prefill-decode, must be > 1, enables batched CPU/GPU prefill");
        out.println("  --print-taskgraph-chain       print every TaskGraph of the GPU plan and when it runs (stderr)");
        out.println();
   */

  /**
   * Ask this question of the LLM.
   */
  @Option(names = {"-p", "--prompt"}, defaultValue = "Please count the number of atoms of each element in the system.",
      description = "Ask this question to the LLM.")
  private String prompt = null;

  /**
   * Include the molecular atom list in the prompt.
   */
  @Option(names = {"-i", "--includeMolecule"}, defaultValue = "false",
      description = "Include the molecular atom list in the prompt.")
  private boolean includeMolecule = false;

  /**
   * Stream tokens as they are generated. This may cause encoding artifacts for non-ASCII text.
   */
  @Option(names = {"--stream"}, defaultValue = "false",
      description = "Stream tokens as they are generated. This may cause encoding artifacts for non-ASCII text.")
  private boolean stream = true;

  /**
   * print ALL tokens to stderr, if true, recommended to set --stream=false, default false
   */
  @Option(names = {"--echo"}, defaultValue = "false",
      description = "Print ALL tokens to stderr, if true, recommended to set --stream=false, default false.")
  private boolean echo = false;

  /**
   * Sampling temperature. The Llama Instruct default is 0.3.
   */
  @Option(names = {"--temperature"}, defaultValue = "0.3",
      description = "Sampling temperature (0 for greedy decoding).")
  private double temperature = 0.3;

  /**
   * Nucleus sampling probability. The Llama Instruct default is 0.95.
   */
  @Option(names = {"--top-p"}, defaultValue = "0.95",
      description = "Nucleus sampling probability.")
  private double topP = 0.95;

  /**
   * Use TornadoVM to run the model.
   *
   * Note that "-Djitllm.kvcache.fp32=true" is needed on MacOS.
   */
  @Option(names = {"-t", "--use-tornadovm"}, defaultValue = "false",
      description = "Use TornadoVM to run the model.")
  private boolean useTornadoVM = false;

  /**
   * Path to a GGUF model file.
   */
  @Option(names = {"-m", "--model"}, required = true,
      description = "Path to a GGUF model file.")
  private String ggufModel = null;

  /**
   * Context capacity, including prompt and generated tokens.
   * If the prompt alone exceeds this limit, no tokens will be generated.
   */
  @Option(names = {"-c", "--ctx-size"}, defaultValue = "1024",
      description = "Context capacity, including prompt and generated tokens. If the prompt alone exceeds this limit, no tokens will be generated.")
  private int contextSize = 1024;

  /**
   * The final argument is a PDB or XYZ coordinate file.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = "An atomic coordinate file in PDB or XYZ format.")
  private String filename = null;

  public Chat() {
    super();
  }

  public Chat(FFXBinding binding) {
    super(binding);
  }

  public Chat(String[] args) {
    super(args);
  }

  @Override
  public Chat run() {
    // Init the context and bind variables.
    if (!init()) {
      return this;
    }

    // Load the MolecularAssembly.
    activeAssembly = getActiveAssembly(filename);
    if (activeAssembly == null) {
      logger.info(helpString());
      return this;
    }

    // Set the filename.
    filename = activeAssembly.getFile().getAbsolutePath();

    Atom[] atoms = activeAssembly.getAtomArray();

    String promptWithMolecule = prompt;
    if (includeMolecule) {
      promptWithMolecule = promptWithMolecule(prompt, atoms, includeMolecule);
      logger.info(promptWithMolecule);
    }

    String[] args = new String[]{
        "run",
        "--model", ggufModel,
        "--instruct",
        "--prompt", promptWithMolecule,
        "--stream", Boolean.toString(stream),
        "--echo", Boolean.toString(echo),
        "--temperature", Double.toString(temperature),
        "--top-p", Double.toString(topP),
        "--ctx-size", Integer.toString(contextSize),
        "--use-tornadovm", Boolean.toString(useTornadoVM),
    };

    Options options = Options.parseOptions(args);
    long startedNs = System.nanoTime();
    ModelOptions modelOptions =
        new ModelRunConfig(
            options.modelPath(),
            options.contextLength(),
            options.useTornadovm())
            .modelOptions();

    try (LocalModel model = LocalModels.load(options.modelPath(), modelOptions)) {
      long modelLoadNs = System.nanoTime() - startedNs;
      guardDeviceSample(model, options);
      try (GenerationSession session = ((TextGenerationModel) model).newSession()) {
        if (StartupDiagnostics.verbose()) {
          String sampling =
              options.temperature() == 0
                  ? "greedy"
                  : String.format(
                  Locale.ROOT,
                  "temperature %.3f / top-p %.3f / seed %d",
                  options.temperature(),
                  options.topp(),
                  options.seed());
          logger.info(StartupDiagnostics.render(
                  model,
                  session.prepare(),
                  sampling,
                  modelOptions,
                  modelLoadNs,
                  startedNs));
        }
        if (options.interactive()) {
          runInteractive(session, options);
        } else {
          runSingleInstruction(session, options);
        }
      }
    } catch (Exception e) {
      logger.log(Level.SEVERE, "Error running Chat command", e);
    }

    return this;
  }

  /**
   * Adds the molecular atom list to a prompt when requested.
   *
   * @param prompt user-provided prompt
   * @param atoms atoms in the active molecular assembly
   * @param includeMolecule whether to append the atom list
   * @return the prompt to send to the LLM
   */
  static String promptWithMolecule(String prompt, Atom[] atoms, boolean includeMolecule) {
    if (!includeMolecule) {
      return prompt;
    }
    StringBuilder moleculePrompt = new StringBuilder(prompt)
        .append(System.lineSeparator())
        .append(System.lineSeparator())
        .append("Atoms in the system:")
        .append(System.lineSeparator());
    for (Atom atom : atoms) {
      moleculePrompt.append(atom).append(System.lineSeparator());
    }
    return moleculePrompt.toString();
  }

  /**
   * On-device greedy sampling ({@code -Djitllm.deviceSample=true}) keeps the logits on the GPU
   * and returns only the argmax token id. It is only valid on the GPU FP16 greedy path for the
   * models whose decode loop reads {@code state.workspace.sampledToken} (Llama / Mistral /
   * Qwen3). For any other configuration the host still needs the full logits row, so the flag is
   * cleared here.
   *
   * <p>Must run after the model is loaded and before a session is opened: the property is read
   * when the session builds its execution plan.
   */
  private static void guardDeviceSample(LocalModel model, Options options) {
    if (!Boolean.getBoolean("jitllm.deviceSample")) {
      return;
    }
    boolean greedy = options.temperature() == 0.0f;
    boolean fp16 = "FP16".equals(model.info().computeType().name());
    String architecture = model.info().architecture();
    boolean wiredLoop =
        architecture.equals("llama")
            || architecture.equals("mistral")
            || architecture.equals("qwen3");
    if (!(options.useTornadovm() && greedy && fp16 && wiredLoop)) {
      logger.warning(
          "[deviceSample] ignored — requires GPU + greedy (temperature 0) + FP16 + Llama/Mistral/Qwen3");
      System.clearProperty("jitllm.deviceSample");
    }
  }

  /** The request shape both modes share; only the prompt and system prompt differ per turn. */
  private static GenerationRequest.Builder request(Options options) {
    return GenerationRequest.builder()
        .maxNewTokens(options.maxNewTokens())
        .temperature(options.temperature())
        .topP(options.topp())
        .seed(options.seed());
  }

  private static void runSingleInstruction(GenerationSession session, Options options) {

    logger.fine(" Running single instruction...");
    logger.fine(" Prompt: " + options.prompt());
    logger.fine(" System prompt: " + options.systemPrompt());
    logger.fine(" Stream: " + options.stream());
    logger.fine(" Echo: " + options.echo());
    logger.fine(" Max tokens: " + options.maxTokens());
    logger.fine(" Use TornadoVM: " + options.useTornadovm());

    GenerationRequest.Builder builder =
        request(options).prompt(options.prompt()).systemPrompt(options.systemPrompt());
    if (options.stream()) {
      builder.onEvent(event -> {
        logger.info(event.text());
      });
    }
    GenerationResult result = session.generate(builder.build());
    if (options.stream()) {
      logger.info(" Streaming complete.");
    }
    if (!options.stream()) {
      logger.info(" Generation complete.");
      logger.info(result.text());
    }
    if (result.finishReason() == FinishReason.CONTEXT_FULL && result.generatedTokens() == 0) {
      // The prompt alone filled the capacity --max-tokens sized: nothing was generated,
      // and printing a zero-token metrics block was the only sign of it.
      throw new IllegalArgumentException(
          contextFullMessage(result.promptTokens(), options.maxTokens()));
    }
    if (SHOW_PERF_INTERACTIVE) {
      RunMetrics.printMetrics();
    }
  }

  /**
   * The diagnostic for a prompt that leaves no room to generate: {@code --max-tokens} is the
   * capacity in positions this run was sized for, prompt included.
   */
  static String contextFullMessage(int promptTokens, int maxTokens) {
    return "the prompt is "
        + promptTokens
        + " tokens and --ctx-size "
        + maxTokens
        + " is the whole capacity (prompt plus generated tokens), so nothing could be"
        + " generated; pass --ctx-size larger than the prompt";
  }

  /**
   * The chat loop. The session carries the conversation, so each turn sends only the new user
   * text; the system prompt goes with the first turn and is retained from there.
   */
  private static void runInteractive(GenerationSession session, Options options) {
    Scanner in = new Scanner(System.in);
    boolean firstTurn = true;
    while (true) {
      logger.info("> ");
      if (!in.hasNextLine()) {
        break;
      }
      String userText = in.nextLine();
      if (userText.equals("quit") || userText.equals("exit")) {
        break;
      }

      GenerationRequest.Builder builder = request(options).prompt(userText);
      if (firstTurn) {
        builder.systemPrompt(options.systemPrompt());
        firstTurn = false;
      }
      if (options.stream()) {
        builder.onEvent(event -> {
          logger.info(event.text());
        });
      }

      GenerationResult result = session.generate(builder.build());
      if (options.stream()) {
        logger.info(""); // newline after the streamed output
      } else {
        logger.info(result.toString());
      }

      if (result.finishReason() == FinishReason.CONTEXT_FULL) {
        logger.warning("\n"
                + contextFullMessage(result.promptTokens(), options.maxTokens())
                + " (or start a new session)");
        break;
      }
      if (SHOW_PERF_INTERACTIVE) {
        RunMetrics.printMetrics();
      }
    }
  }

}
