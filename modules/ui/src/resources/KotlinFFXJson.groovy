import groovy.io.FileType

// Constant header of the Kotlin ffx.json file.
header = """{
  "properties": {
    "v":"1.0.0-beta"
  },
  "link": "https://ffx.biochem.uiowa.edu",
  "repositories": [
    "file:/ffx/lib"
  ],
  "imports": [
    "ffx.Main",
    "ffx.Main.*"
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
dir.eachFileRecurse (FileType.FILES) { file ->
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
File json = new File("binder/ffx.json")
json.write(sb.toString())

