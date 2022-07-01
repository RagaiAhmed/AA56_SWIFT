import XCTest
@testable import AA56

final class temp_swiftTests: XCTestCase {
    func testExample() throws {
        print("fuf");
        var res = calcPolar(1656626400,0);

	print("Sun on 1/7/2022 at 00:00")
	print("Longitude ",res.lon)
	print("Declination ",res.dec)
	print("Distance ",res.r)
    }
}
