// Miscellaneous Utils

def logWithColor(message, color) {
    log.info "${color}${message}\u001B[0m"
}
def logSuccess(message) {
    logWithColor(message, "\u001B[32m")  // Green
}
def logWarning(message) {
    logWithColor(message, "\u001B[33m")  // Yellow
}
def logError(message) {
    logWithColor(message, "\u001B[31m")  // Red
}
def logInfo(message) {
    logWithColor(message, "\u001B[36m")  // Cyan
}