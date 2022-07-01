import XCTest
@testable import AA56


final class temp_swiftTests: XCTestCase {
    func testExample() throws {
	
	var path = UnsafeMutablePointer<CChar>(mutating: "aa.ini".utf8String)
	var status = initCalc(path)
	print("Init Status (0 means successful) : ",status)
        var res = calcPolar(1656626400,0)

	print("Sun on 1/7/2022 at 00:00")
	print("Longitude ",res.lon)
	print("Declination ",res.dec)
	print("Distance ",res.r)
        
	status = initCalc(path)
	print("Init Status (0 means successful) : ",status)
        res = calcPolar(1656626400,1);

	print("Mercury on 1/7/2022 at 00:00")
	print("Longitude ",res.lon)
	print("Declination ",res.dec)
	print("Distance ",res.r)
	
	status = initCalc(path)
	print("Init Status (0 means successful) : ",status)
        res = calcPolar(1656626400,2);

	print("Venus on 1/7/2022 at 00:00")
	print("Longitude ",res.lon)
	print("Declination ",res.dec)
	print("Distance ",res.r)
	
		
	status = initCalc(path)
	print("Init Status (0 means successful) : ",status)
        res = calcPolar(1656626400,3);

	print("Moon on 1/7/2022 at 00:00")
	print("Longitude ",res.lon)
	print("Declination ",res.dec)
	print("Distance ",res.r)
	
		
	status = initCalc(path)
	print("Init Status (0 means successful) : ",status)
        res = calcPolar(1656626400,4);

	print("Mars on 1/7/2022 at 00:00")
	print("Longitude ",res.lon)
	print("Declination ",res.dec)
	print("Distance ",res.r)
	
		
	status = initCalc(path)
	print("Init Status (0 means successful) : ",status)
        res = calcPolar(1656626400,5);

	print("Jupiter on 1/7/2022 at 00:00")
	print("Longitude ",res.lon)
	print("Declination ",res.dec)
	print("Distance ",res.r)
	
		
	status = initCalc(path)
	print("Init Status (0 means successful) : ",status)
        res = calcPolar(1656626400,6);

	print("Saturn on 1/7/2022 at 00:00")
	print("Longitude ",res.lon)
	print("Declination ",res.dec)
	print("Distance ",res.r)
	
		
	status = initCalc(path)
	print("Init Status (0 means successful) : ",status)
        res = calcPolar(1656626400,7);

	print("Uranus on 1/7/2022 at 00:00")
	print("Longitude ",res.lon)
	print("Declination ",res.dec)
	print("Distance ",res.r)
	
		
	status = initCalc(path)
	print("Init Status (0 means successful) : ",status)
        res = calcPolar(1656626400,8);

	print("Neptune on 1/7/2022 at 00:00")
	print("Longitude ",res.lon)
	print("Declination ",res.dec)
	print("Distance ",res.r)
	
	
		
	status = initCalc(path)
	print("Init Status (0 means successful) : ",status)
        res = calcPolar(1656626400,9);

	print("Pluto on 1/7/2022 at 00:00")
	print("Longitude ",res.lon)
	print("Declination ",res.dec)
	print("Distance ",res.r)
    }
}
