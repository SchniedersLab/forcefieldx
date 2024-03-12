// Generate the Kotlin ffx.json file for the Jupyter notebook.

// Constant header of the Kotlin ffx.json file.
header = """{
  "properties": {
    "v":"1.0.0-beta"
  },
  "link": "https://ffx.biochem.uiowa.edu",
  "repositories": [
    "file:../lib"
  ],
  "imports": [
    "ffx.Main",
    "ffx.Main.*",
    "ffx.utilities.DownloadUtils.*"
  ],
  "init": [
   "System.setProperty(\\\"java.awt.headless\\\",\\\"true\\\")",
   "System.setProperty(\\\"j3d.rend\\\",\\\"noop\\\")"
  ],
  "dependencies": [
"""
StringBuilder sb = new StringBuilder(header)

// Collect the classpath entries from the lib directory.
def dir = new File("lib")
dir.listFiles().sort{ it.name }.each { file ->
  sb.append("    \"" + file.getName() + "\",\n")
}
last = sb.lastIndexOf(",")
sb.deleteCharAt(last)

// Append the constant footer.
footer = """  ]
}
"""
sb.append(footer)

// Write the Kotlin ffx.json file.
File json = new File("ipynb-kotlin/ffx.json")
json.write(sb.toString())
