// Logging Utils

reset = "\u001B[0m"
red = "\u001B[31m"
green = "\u001B[32m"
yellow = "\u001B[33m"
cyan = "\u001B[36m"

def logWithColor(message, color) {
    log.info "${color}${message}${reset}"
}
def logSuccess(message) {
    logWithColor(message, green)
}
def logWarning(message) {
    logWithColor(message, yellow)
}
def logError(message) {
    logWithColor(message, red)
}
def logInfo(message) {
    logWithColor(message, cyan)
}
def coloredTitle(sep = "") {
    return [
        "\uD83E\uDDFC ", // soap emoji
        "${yellow}F${reset}",
        "${yellow}F${reset}",
        "${yellow}P${reset}",
        "${green}E${reset}",
        "${cyan}R${reset}",
        "${cyan}A${reset}",
        "${cyan}S${reset}",
        "${cyan}E${reset}",
    ].join(sep)
}

// File Utils

def getAbsolute(path) {
    return new File(path).absolutePath
}

def mkdirs(path) {
    path = new File(new File(path).absolutePath)
    if (!path.exists()) {
        path.mkdirs()
    }
}

def logDirTree(String dir) {
    def dirFile = new File(dir)
    if (!dirFile.exists()) {
        logWarning("Directory does not exist: ${dir}")
        return
    }
    
    def tree = new StringBuilder()
    tree.append("\nOutput directory structure:\n")
    tree.append(dir + "\n")
    
    createDirTree(dirFile, "", tree)
    logInfo(tree.toString())
}

def createDirTree(File dir, String indent, StringBuilder tree) {
    File[] files = dir.listFiles()
    if (files == null) return
    
    for (int i = 0; i < files.length; i++) {
        File file = files[i]
        boolean isLast = (i == files.length - 1)
        
        tree.append(indent)
        tree.append(isLast ? "└── " : "├── ")
        tree.append(file.getName() + "\n")
        
        if (file.isDirectory()) {
            String newIndent = indent + (isLast ? "    " : "│   ")
            createDirTree(file, newIndent, tree)
        }
    }
}