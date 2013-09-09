import httplib
import json
import sys

def CreateJSON():
    data = []
    data.append({'name': 'glycine', 
                 'InChI':'InChI=1S/C2H5NO2/c3-1-2(4)5/h1,3H2,(H,4,5)',
                 'pKas': [2.31, 9.24]})
    
    data.append({'name': 'ATP', 
                 'InChI':'InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p-3',
                 'pKas': [-7.03, -3.85, 0.87, 2.31, 2.83, 3.38, 5.13, 7.4, 12.6]})
    
    data.append({'name': 'CO2', 
                 'InChI':'InChI=1S/CO2/c2-1-3',
                 'pKas': []})
    
    data.append({'name': '3-Ketoarabinitol', 
                 'InChI':'InChI=1S/C5H10O5/c6-1-3(8)5(10)4(9)2-7/h3-4,6-9H,1-2H2',
                 'pKas': [-8.64, -3.87, -3.87, -3.01, -3.01, 11.29, 13.13, 15.12, 15.86, 18.62]})
    
    return json.dumps(data)

def Test():
    http = httplib.HTTPConnection('132.77.81.102', 80)
    #http = httplib.HTTPConnection('localhost', 8000)

    url = '/inchi'
    headers = {"Content-Type" : "application/json"}
    body = CreateJSON()
    http.request(method='GET', url=url, headers=headers, body=body)
    response = http.getresponse()
    
    sys.stderr.write(response.fp.read())
    
if __name__ == "__main__":
    Test()