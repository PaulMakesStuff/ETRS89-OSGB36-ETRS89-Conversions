const os = require('express').Router();
const bp = require('body-parser');
const fs = require('fs');

os.use(bp.json());

os.post("/lltoen", function (req, res) {

    let lat = req.body.lat * (Math.PI/180); //this.degrees_to_radians(req.body.lat);
    let lon = req.body.lon * (Math.PI/180); //this.degrees_to_radians(req.body.lon);

    // WGS84 constants
    const a = 6378137;
    const b = 6356752.3141;
    const e2 = (a**2-b**2)/a**2;

    // National Grid projection constants
    const F0 = 0.9996012717;
    const lat0 = 49 * (Math.PI/180); //this.degrees_to_radians(49);
    const lon0 = -2 * (Math.PI/180); //this.degrees_to_radians(-2);
    const E0 = 400000;
    const N0 = -100000;

    // calculation steps
    let n = (a-b)/(a+b);
    let v = a*F0*(1-(e2*(Math.sin(lat)**2)))**-0.5;
    let p = a*F0*(1-e2)*(1-(e2*(Math.sin(lat)**2)))**-1.5;
    let n2 = v/p-1;
    let M = b*F0*((1+n+5/4*n**2+5/4*n**3)*(lat-lat0)-(3*n+3*n**2+21/8*n**3)*Math.sin(lat-lat0)*Math.cos(lat+lat0)+(15/8*n**2+15/8*n**3)*Math.sin(2*(lat-lat0))*Math.cos(2*(lat+lat0))-35/24*n**3*Math.sin(3*(lat-lat0))*Math.cos(3*(lat+lat0)));
    let I = M + N0;
    let II = v/2*Math.sin(lat)*Math.cos(lat);
    let III = v/24*Math.sin(lat)*Math.cos(lat)**3*(5-Math.tan(lat)**2+9*n2);
    let IIIA = v/720*Math.sin(lat)*Math.cos(lat)**5*(61-58*Math.tan(lat)**2+Math.tan(lat)**4);
    let IV = v * Math.cos(lat);
    let V = v/6*Math.cos(lat)**3*(v/p-Math.tan(lat)**2);
    let VI = v/120*Math.cos(lat)**5*(5-18*Math.tan(lat)**2+Math.tan(lat)**4+14*n2-58*(Math.tan(lat)**2*n2));
    let northings = I+II*(lon-lon0)**2+III*(lon-lon0)**4+IIIA*(lon-lon0)**6;
    let eastings = E0+IV*(lon-lon0)+V*(lon-lon0)**3+VI*(lon-lon0)**5;
    
    // perform OSTM15 shift
    let e_idx = parseInt(eastings/1000);
    let n_idx = parseInt(northings/1000);
    let x0 = e_idx * 1000;
    let y0 = n_idx * 1000;

    let corners = [[e_idx, n_idx], [e_idx+1, n_idx], [e_idx+1, n_idx+1], [e_idx, n_idx+1]];
    let record_nums = [];
    for(let i = 0; i < 4; i++){
        record_nums.push(corners[i][0]+(corners[i][1]*701)+1);
    }

    fs.readFile(__dirname + '/OSTN15.csv', 'utf8', function(err, data) {
        if(err){return console.log(err);}
        const lines = data.split("\n");
        let dx = eastings - x0;
        let dy = northings - y0;
        let t = dx / 1000;
        let u = dy / 1000;
        let sea = [];
        let sna = [];
        for(let i =0; i < 4; i++){
            sea.push(lines[record_nums[i]].split(",")[0]);
            sna.push(lines[record_nums[i]].split(",")[1]);
        }
        let se = (1-t)*(1-u)*sea[0]+(t)*(1-u)*sea[1]+(t)*(u)*sea[2]+(1-t)*(u)*sea[3];
        let sn = (1-t)*(1-u)*sna[0]+(t)*(1-u)*sna[1]+(t)*(u)*sna[2]+(1-t)*(u)*sna[3];
        res.json({e: eastings+se, n:northings+sn});
    });       
});


module.exports = os;
