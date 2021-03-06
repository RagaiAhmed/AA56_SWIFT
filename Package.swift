// swift-tools-version:5.1
// The swift-tools-version declares the minimum version of Swift required to build this package.
import PackageDescription

let package = Package(
    name: "AA56",
    platforms: [
        .macOS(.v10_11), .iOS(.v11),
    ],
    products: [
        // Products define the executables and libraries produced by a package, and make them visible to other packages.
        .library(
            name: "AA56",
            type: .dynamic,
            targets: ["AA56"]
        )
    ],
    targets: [
        // Targets are the basic building blocks of a package. A target can define a module or a test suite.
        // Targets can depend on other targets in this package, and on products in packages which this package depends on.
        .target(
            name: "AA56",
            dependencies: [],
            path: "Sources/AA56"
        ),
        .testTarget(
            name: "temp_swiftTests",
            dependencies: ["AA56"])
    ]
)
