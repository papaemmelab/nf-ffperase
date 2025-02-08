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