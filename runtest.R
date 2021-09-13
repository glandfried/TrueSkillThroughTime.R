if (!require("RUnit", quietly = TRUE)) {
  stop("Package Runit is not found.") 
}

# Path to the unit tests folder in the package
dir <- system.file("./R", package="trueskillthroughtime")


# Define RUnit test suite
suite <- defineTestSuite(name=paste("trueskillthroughtime", "RUnit Tests"),
                         dirs="unitTests",
                         testFileRegexp = "runit.R",
                         testFuncRegexp = "^test_+",
                         rngKind="default",
                         rngNormalKind="default")
                         
# Run tests
result <- runTestSuite(suite)

# Display result tests on the console
printTextProtocol(result)
