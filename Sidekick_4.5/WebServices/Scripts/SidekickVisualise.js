function maximum (myArray) {
var max = myArray[0];
var len = myArray.length;
for (var i = 1; i < len; i++) if (myArray[i] > max) max = myArray[i];
return max;
}

function maximum2d (myArray) {
	var max = myArray[0][0];
	var len = myArray.length;
	var rlen = myArray[0].length;
	for (var i = 0; i < len; i++) {
		for (var j = 0; j < len; j++) {
			if (myArray[i][j] > max) max = myArray[i][j];
			}
		}
	return max;
}

function lerp(min,max,value) {
//gives the fractional distance the value is between the max and min
	value = value - min;
	max = max - min;
	return value/max;
}

function colormap(svalue,min,max) {
	value = lerp(min,max,svalue);
	//return "rgb(0,0,0)";
	//0.875
	if (value>0.875) {
		blue=0;
		green=0;
		red=Math.round(256-128*lerp(0.875,1,value));
	} else if (value<0.125) {
		red=0;
		green=0;
		blue=Math.round(128+128*lerp(0,0.125,value));
	} else if (value<0.5) {
		red=0;
		green=Math.round(256*lerp(0.125,0.5,value));
		blue=Math.round(256-256*lerp(0.125,0.5,value));
	}else{
		blue=0;
		green=Math.round(256-256*lerp(0.5,0.875,value));
		red=Math.round(256*lerp(0.5,0.875,value));
	}
	//0.5

	//0.125	
	
	
	return "rgb("+red+","+green+","+blue+")";
	}

function drawHistogram(data,bins,whichDiv,maxvalue,cbar,colorby,logbin) {
	//var whichDiv = "canvas";
	//alert(data);
	var canvas = document.getElementById(whichDiv);
	var width = canvas.width;
	var height = canvas.height;

	if (cbar){
		wmod=0.7;
		} else {
		wmod=0.8;
		}
	if (maxvalue==-1){
		maxvalue=maximum(data);
		}
	if (logbin) {
		maxvalue=Math.log(maxvalue);
		}
	if (canvas.getContext) {
		var ctx = canvas.getContext("2d");
		
		//background
		ctx.fillStyle="rgb(256,256,256)";
		ctx.fillRect (0,0,width,height);
		
		ctx.linewidth = 2;
		ctx.lineCap = "square";
		if (width%2 == 1) {
			graphLedgeX = width*0.1;
			graphTedgeY = height*0.1;
		} else {
			graphLedgeX = width*0.1+0.5;
			graphTedgeY = height*0.1+0.5;
			}
		ctx.strokeRect(graphLedgeX,graphTedgeY,width*wmod,height*0.8);
		
		
		//draw bars
		
		bar_width = width*wmod/(data.length);
		
		ctx.save();
		
		ctx.translate(graphLedgeX,graphTedgeY);
		
		ctx.save();
		//normalise in future? a shared max?
		transform = height*0.8/maxvalue;
		//ctx.fillStyle="rgb(0,0,0)";
		for (id in data) {
			cvalue = data[id];
			if (logbin){
				if ( cvalue > 0 ) {
					cvalue = Math.log(cvalue);
					}
				}
			if (colorby=="position"){
				ctx.fillStyle=colormap(bins[id],bins[0], bins[bins.length-2]);
			} else if (colorby=="black") {
				ctx.fillStyle='rgb(0,0,0)';
			} else {
				ctx.fillStyle=colormap(cvalue, 0, maxvalue);
			}
			ctx.strokeRect(0,height*0.8-(transform*cvalue),bar_width,transform*cvalue);			
			ctx.fillRect(0,height*0.8-(transform*cvalue),bar_width,transform*cvalue);
			ctx.translate(bar_width,0);
			
		}
		
		ctx.restore();ctx.save();
		
		//now for the tick marks and legends
		//X
		ctx.translate(0,height*0.8);
		ctx.strokeStyle="rgb(0,0,0)";
		ctx.lineWidth=1;

		ctx.fillStyle="black";
		ctx.textAlign="center";
		fontsize = height/15;
		ctx.font = fontsize + "px Helvetica,sans-serif";

		for (id in bins){
			if (id%5==0){
				ctx.beginPath();
				ctx.moveTo(0,0);
				ctx.lineTo(0,width/80);
				
				ctx.stroke();
				
				ctx.fillText(bins[id],0,width/80+fontsize);
				

				//ctx.arc(0,i*12.5,5,0,Math.PI*2,true);
				//ctx.beginPath();
				//ctx.arc(0,00,5,0,Math.PI*2,true);
				//ctx.stroke();
				}
			ctx.translate(bar_width,0);
			}
		ctx.restore();
		
		ctx.save();

		//Y
		ctx.textAlign="right";
		ctx.fillStyle="black";
		ctx.translate(-width/80,0);
		ctx.translate(0,height*0.8);
		lpu=(height*0.8)/maxvalue;
		//let's work on the assumption that we want about 5
		divisions= maxvalue/5;
		//alert(divisions);
		//we probably want these rounded to the nearest 10^x
		oom=Math.pow(10,Math.round((Math.log(divisions)-0.5)/2.3));
		//alert(oom);
		divisions=Math.round(divisions/oom)*oom;
		//alert(divisions);
		//alert(data.max());
		for (i=0;i< maxvalue ;i=i+divisions){
				ctx.beginPath();
				ctx.moveTo(0,0);
				ctx.lineTo(width/80,0);				
				ctx.stroke();
				ctx.fillText(i/oom,-fontsize/2,fontsize/3);
				//alert(i)
				ctx.translate(0,-divisions*lpu);
			}
		//ctx.fillText("y=y/"+oom,fontsize*3,fontsize/3);

		ctx.restore();
		ctx.save();
		//y key
		ctx.fillStyle="black";
		ctx.textAlign="left";
		fontsize = height/15;
		ctx.font = fontsize + "px Helvetica,sans-serif";
		if (logbin) {
			if (oom == 1) {
				ctx.fillText("y=ln(y)",-fontsize*2,-fontsize/2);
				} else {
				ctx.fillText("y=ln(y)/"+oom,-fontsize*2,-fontsize/2);
				}
			} else {
			if (oom == 1) {
				} else {
				ctx.fillText("y=y/"+oom,-fontsize*2,-fontsize/2);
				}
			}
		ctx.restore();
		//colorbar
		if (cbar){
		
		ctx.translate(width*0.75,0);
		ctx.strokeRect(0,0,width*0.05,height*0.8);

		lingrad=ctx.createLinearGradient(0,0,0,height*0.8);
		lingrad.addColorStop(0,'rgb(128,0,0)');
		lingrad.addColorStop(0.125,'rgb(256,0,0)');
		lingrad.addColorStop(0.5,'rgb(0,256,0)');
		lingrad.addColorStop(0.875,'rgb(0,0,256)');
		lingrad.addColorStop(1,'rgb(0,0,128)');
		ctx.fillStyle=lingrad;
		ctx.fillRect(0,0,width*0.05,height*0.8);
		}

	ctx.restore();
	}
}

function draw2dHistogram(data,xbins,ybins,whichDiv,maxvalue,cbar,colorby,logbin) {
	//var whichDiv = "canvas";
	//alert(data);
	var canvas = document.getElementById(whichDiv);
	var width = canvas.width;
	var height = canvas.height;
	//alert(data[0])
	if (cbar){
		wmod=0.7;
		} else {
		wmod=0.8;
		}
	if (maxvalue==-1){
		maxvalue=maximum2d(data);
		}
	if (logbin) {
		maxvalue = Math.log(maxvalue);
		}
	if (canvas.getContext) {
		var ctx = canvas.getContext("2d");
		
		//background
		ctx.fillStyle="rgb(256,256,256)";
		ctx.fillRect (0,0,width,height);
		
		ctx.linewidth = 2;
		ctx.lineCap = "square";
		if (width%2 == 1) {
			graphLedgeX = width*0.1;
			graphTedgeY = height*0.1;
		} else {
			graphLedgeX = width*0.1+0.5;
			graphTedgeY = height*0.1+0.5;
			}
		ctx.strokeRect(graphLedgeX,graphTedgeY,width*wmod,height*0.8);
		
		
		//draw squares
		
		bar_width = width*wmod/(data.length);
		bar_height = height*0.8/(data[0].length);
		
		ctx.save();
		
		ctx.translate(graphLedgeX,graphTedgeY);
		
		ctx.save();
		//normalise in future? a shared max?
		transform = height*0.8/maxvalue;
		ctx.translate(width*wmod,height*0.8);
		ctx.rotate(Math.PI);
		ctx.translate(width*wmod,0);
		//ctx.fillStyle="rgb(0,0,0)";
		for (id in data) {
			ctx.save()
			for (jd in data[id]) {
				if (data[id][jd] > 0) {
					if (logbin) {
						ctx.fillStyle=colormap(Math.log(data[id][jd]), 0, maxvalue);
						ctx.strokeStyle=colormap(Math.log(data[id][jd]), 0, maxvalue);
						} else {
						ctx.fillStyle=colormap(data[id][jd], 0, maxvalue);
						ctx.strokeStyle=colormap(data[id][jd], 0, maxvalue);
						}
					} else {
					ctx.fillStyle="rgb(256,256,256)";
					ctx.strokeStyle="rgb(256,256,256)";	
					}			
				ctx.fillRect(-bar_width,0,bar_width,bar_height);
				ctx.strokeRect(-bar_width,0,bar_width,bar_height);
				ctx.translate(0,bar_height);
			}
			ctx.restore()
			ctx.translate(-bar_width,0);
		}
		
		ctx.restore();
		
		ctx.save();
		
		//Redraw histogram edges
		ctx.strokeRect(0,0,width*wmod,height*0.8);
		
		//X
		ctx.translate(0,height*0.8);
		ctx.strokeStyle="rgb(0,0,0)";
		ctx.lineWidth=1;

		ctx.fillStyle="black";
		ctx.textAlign="center";
		fontsize = height/15;
		ctx.font = fontsize + "px Helvetica,sans-serif";

		for (id in xbins){
			if (id%5==0){
				ctx.beginPath();
				ctx.moveTo(0,0);
				ctx.lineTo(0,width/80);
				
				ctx.stroke();
				
				ctx.fillText(xbins[id],0,width/80+fontsize);
				

				//ctx.arc(0,i*12.5,5,0,Math.PI*2,true);
				//ctx.beginPath();
				//ctx.arc(0,00,5,0,Math.PI*2,true);
				//ctx.stroke();
				}
			ctx.translate(bar_width,0);
			}
		ctx.restore();
		
		//Y
		ctx.translate(0,height*0.8);
		ctx.strokeStyle="rgb(0,0,0)";
		ctx.lineWidth=1;

		ctx.fillStyle="black";
		ctx.textAlign="right";
		fontsize = height/15;
		ctx.font = fontsize + "px Helvetica,sans-serif";

		for (id in ybins){
			if (id%5==0){
				ctx.beginPath();
				ctx.moveTo(0,0);
				ctx.lineTo(-width/80,0);
				
				ctx.stroke();
				
				ctx.fillText(ybins[id],-width/80,fontsize/3);
				

				//ctx.arc(0,i*12.5,5,0,Math.PI*2,true);
				//ctx.beginPath();
				//ctx.arc(0,00,5,0,Math.PI*2,true);
				//ctx.stroke();
				}
			ctx.translate(0,-bar_height);
			}
		
		ctx.restore();
				
	}
}


function helix_cartoon (whichDiv, position, tilt, rotation, sidechains, sequence, rotateView,thickness) {
	var canvas = document.getElementById(whichDiv);
	var width = canvas.width;
	var height = canvas.height;
	sequence_length = sequence.length;

	var amino_color_lookup = new Array();
	amino_color_lookup["L"] = "white";
	amino_color_lookup["I"] = "white";
	amino_color_lookup["V"] = "white";
	amino_color_lookup["A"] = "white";
	amino_color_lookup["M"] = "white";
	amino_color_lookup["P"] = "white";
	amino_color_lookup["G"] = "white";
	amino_color_lookup["S"] = "rgb(60,155,210)";
	amino_color_lookup["T"] = "rgb(60,155,210)";
	amino_color_lookup["N"] = "rgb(60,155,210)";
	amino_color_lookup["Q"] = "rgb(60,155,210)";
	amino_color_lookup["C"] = "rgb(60,155,210)";
	amino_color_lookup["W"] = "green";
	amino_color_lookup["Y"] = "green";
	amino_color_lookup["F"] = "green";
	amino_color_lookup["R"] = "blue";
	amino_color_lookup["K"] = "blue";
	amino_color_lookup["H"] = "blue";
	amino_color_lookup["D"] = "red";
	amino_color_lookup["E"] = "red";
	
	if (rotateView) {
		tilt = -tilt;
		//alert(rotation);
		rotation = (180 + rotation)%360;
		//alert(rotation);
		}
	
	if (canvas.getContext) {
	
		
		var ctx = canvas.getContext("2d");
		ctx.save()
		
		//background
		ctx.fillStyle="rgb(60,155,210)";
		ctx.fillRect(0,0,width,height);
		
		//membrane
		//Assume DPPC for now- 41A should be the middle third;
		pixels_per_angstrom = (height/3)/41;
		bilayer_width = thickness * pixels_per_angstrom;
		ctx.fillStyle="rgb(256,256,256)";
		ctx.fillRect(0,(height/2-bilayer_width/2),width,bilayer_width);
		fontsize = height/30;
		ctx.font = fontsize + "px Helvetica,sans-serif";
		
		//membrane center
		//ctx.beginPath();
		//ctx.moveTo(0,height/2);
		//ctx.lineTo(width,height/2);
		//ctx.stroke();
		
		//draw helix
		helix_length = sequence_length*1.5*pixels_per_angstrom;
		helix_radius = 6*pixels_per_angstrom;
		ctx.translate(width/2,height/2);
		pixel_position = position*pixels_per_angstrom;

		//alert([bilayer_width,pixel_position,helix_radius]);
		ctx.translate(0,-pixel_position);
		ctx.rotate(tilt*Math.PI/180);
		
		if (sidechains) {
			//define sidechains as behind or in front of the helix
			//alert(1)
			angle_0 = (360 + rotation - 200)%360;
			angles = [];
			for (var i = 0; i < sequence_length; i++) {
				angles[i] = (angle_0 + i*100)%360;
				}
			//draw sidechains behind helix
			ctx.fillStyle = "white";
			ctx.strokeStyle = "black";
			//alert(angles)
			//alert(2)
			for (id in angles) {
				if ( angles[id] < 180) {
					//ctx.fillstyle = "rgb("+Math.floor(id/(helix_length-1)*256)+",0,0)";
					//if (angles[id] < 180) { ctx.fillStyle = "black"; } else { ctx.fillStyle = "white"; }
					//alert(3)
					ctx.fillStyle = amino_color_lookup[sequence.charAt(id)];
					y_projection = -(-helix_length/2 - 0.75*pixels_per_angstrom + helix_length) + id*1.5*pixels_per_angstrom;
					x_projection = -Math.cos(angles[id]*Math.PI/180)*(helix_radius+2.4*pixels_per_angstrom);
					//alert(3.5)
					ctx.beginPath();
					//alert([x_projection,y_projection,2.4*pixels_per_angstrom,0,2*Math.PI]);
					ctx.arc(x_projection,y_projection,2.4*pixels_per_angstrom,0,2*Math.PI,false);
					ctx.fill();
					//alert(4)
					ctx.beginPath();
					ctx.arc(x_projection,y_projection,2.4*pixels_per_angstrom,0,2*Math.PI,false);
					ctx.stroke();
					ctx.textAlign = "center";
					ctx.fillStyle = "black";
					ctx.fillText((parseInt(id)+1),x_projection,y_projection+fontsize/3);
					//ctx.fillText(sequence.charAt(id),x_projection,y_projection);
					//alert(sequence.charAt(id))
					ctx.fillStyle = "white";
					//alert(5)
					//alert([id,angles[id]]);
					}
				}
			}
		
		//ctx.fillStyle="rgb(100,256,100)";
		lingrad = ctx.createLinearGradient(-helix_radius,0,helix_radius,0);
		if (rotateView) {
			lingrad.addColorStop(0,'rgb(100,256,100)');
			lingrad.addColorStop(0.66,'rgb(100,256,100)');
			lingrad.addColorStop(1,'rgb(256,256,256)');
			} else {
			lingrad.addColorStop(0,'rgb(256,256,256)');
			lingrad.addColorStop(0.33,'rgb(100,256,100)');
			lingrad.addColorStop(1,'rgb(100,256,100)');
			}
		
		ctx.fillStyle = lingrad;
		
		ctx.fillRect(-helix_radius,-helix_length/2,2*helix_radius,helix_length);
		ctx.strokeRect(-helix_radius,-helix_length/2,2*helix_radius,helix_length);
		
		ctx.fillStyle="rgb(100,256,100)";
		//ctx.beginPath();
		//ctx.moveTo(0,-helix_length/2);
		//ctx.arc(0,-helix_length/2,helix_radius,0,2*Math.PI);
		//ctx.fill();
		
		ctx.fillStyle = lingrad;
		//ctx.beginPath();
		//ctx.moveTo(0,-helix_length/2);
		//ctx.arc(0,-helix_length/2,helix_radius,0,2*Math.PI);
		//ctx.stroke();
		
		//ctx.beginPath();
		//ctx.moveTo(0,-helix_length/2);
		//ctx.arc(0,helix_length/2,helix_radius,0,2*Math.PI);
		//ctx.fill();
		
		//ctx.beginPath();
		//ctx.moveTo(0,-helix_length/2);
		//ctx.arc(0,helix_length/2,helix_radius,0,Math.PI);
		//ctx.stroke();

		if (sidechains) {
			//draw sidechains in front of helix
			ctx.fillStyle = "white";
			ctx.strokeStyle = "black";
			//alert(angles)
			for (id in angles) {
				if (angles[id] >= 180) {
					//ctx.fillstyle = "rgb("+Math.floor(id/(helix_length-1)*256)+",0,0)";
					//if (angles[id] < 180) { ctx.fillStyle = "black"; } else { ctx.fillStyle = "white"; }
					ctx.fillStyle = amino_color_lookup[sequence.charAt(id)];
					y_projection = -(-helix_length/2 - 0.75*pixels_per_angstrom + helix_length) + id*1.5*pixels_per_angstrom;
					x_projection = -Math.cos(angles[id]*Math.PI/180)*(helix_radius+2.4*pixels_per_angstrom);
					ctx.beginPath();
					ctx.arc(x_projection,y_projection,2.4*pixels_per_angstrom,0,2*Math.PI,false);
					ctx.fill();
					ctx.beginPath();
					ctx.arc(x_projection,y_projection,2.4*pixels_per_angstrom,0,2*Math.PI,false);
					ctx.stroke();
					//alert([id,angles[id]]);
					ctx.textAlign = "center";
					ctx.fillStyle = "black";
					ctx.fillText((parseInt(id)+1),x_projection,y_projection+fontsize/3);
					ctx.fillStyle = "white";
					}
				}
			}
		
		ctx.restore();
		
		}
	
	}

function read_array ( somestring ) {
	number_array = somestring.split(',');
	for (number_id in number_array) {
		number_array[number_id] = parseInt(number_array[number_id]);
		}
	return number_array;
	}
	
function read_2d_array (somestring) {
	two_d_array = somestring.split('/');
	for (id in two_d_array) {
		two_d_array[id] = read_array(two_d_array[id]);
		}
	return two_d_array;
	}

function read_3d_array (somestring) {
	three_d_array = somestring.split(":");
	for (id in three_d_array) {
		three_d_array[id] = read_2d_array(three_d_array[id]);
		}
	//alert([three_d_array.length,three_d_array[0].length,three_d_array[0][0].length]);
	return three_d_array;
	}


function get_modal_value (histogram, bins) {
	modal_id = 0;
	modal_value = histogram[0];
	for (id in histogram) {
		if (histogram[id]>modal_value) {
			modal_id = id;
			modal_value = histogram[id];
			}
		}
	bin_midpoint = (bins[0]+bins[1])/2-bins[0];
	modal_centerpoint = bins[modal_id] + bin_midpoint;
	//alert(modal_centerpoint);
	return [modal_centerpoint, modal_id];
	}

function get_3d_modal_value (histogram, xbins,ybins,zbins) {
	modal_id = [0,0,0];
	modal_value = histogram[0][0][0];
	for (xid in histogram) {
		for (yid in histogram[xid]) {
			for (zid in histogram[xid][yid]) {			
				if (histogram[xid][yid][zid]>modal_value) {
					modal_id = [xid,yid,zid];
					modal_value = histogram[xid][yid][zid];
					//alert(modal_id);
					//alert(modal_value);
					}
				}
			}
		}
	xbin_midpoint = (xbins[0]+xbins[1])/2-xbins[0];
	ybin_midpoint = (ybins[0]+ybins[1])/2-ybins[0];
	zbin_midpoint = (zbins[0]+zbins[1])/2-zbins[0];
	xmodal_centerpoint = xbins[modal_id[0]] + xbin_midpoint;
	ymodal_centerpoint = ybins[modal_id[1]] + ybin_midpoint;
	zmodal_centerpoint = zbins[modal_id[2]] + zbin_midpoint;
	//alert(modal_centerpoint);
	//alert([xmodal_centerpoint,ymodal_centerpoint,zmodal_centerpoint]);
	//alert(modal_id);
	//alert(modal_value);
	return [xmodal_centerpoint,ymodal_centerpoint,zmodal_centerpoint];
	}

function get_2d_modal_value (histogram,xbins,ybins){
	modal_id = [0,0,0];
	modal_value = histogram[0][0];
	for (xid in histogram) {
		for (yid in histogram[xid]) {
				if (histogram[xid][yid]>modal_value) {
					modal_id = [xid,yid];
					modal_value = histogram[xid][yid];
					//alert(modal_id);
					//alert(modal_value);
					}
			}
		}
	xbin_midpoint = (xbins[0]+xbins[1])/2-xbins[0];
	ybin_midpoint = (ybins[0]+ybins[1])/2-ybins[0];
	xmodal_centerpoint = xbins[modal_id[0]] + xbin_midpoint;
	ymodal_centerpoint = ybins[modal_id[1]] + ybin_midpoint;
	//alert(modal_centerpoint);
	//alert([xmodal_centerpoint,ymodal_centerpoint,zmodal_centerpoint]);
	//alert(modal_id);
	//alert(modal_value);
	return [xmodal_centerpoint,ymodal_centerpoint, modal_id];
	} 

function detailStartup(){
	//document.getElementById('canvas').width=400;
	//document.getElementById('canvas').height=300;
	var data = new Array();
	try {	
	xmlglial=new XMLHttpRequest();
	xmlglial.open("GET",document.getElementById("datasetSelect").value+"?",false);
	xmlglial.send();
	glialDoc=xmlglial.responseXML;
	} catch (e) {
	alert("an error was encountered:"+e) }
	data['sequence'] = glialDoc.getElementsByTagName("sequence")[0].childNodes[0].nodeValue
	data['thickness'] = parseInt(glialDoc.getElementsByTagName("thickness")[0].childNodes[0].nodeValue)
	
	data['position'] = read_array(glialDoc.getElementsByTagName("position")[0].childNodes[0].nodeValue);
	data['position_bin'] = read_array(glialDoc.getElementsByTagName("position_bin")[0].childNodes[0].nodeValue);
	data['tilt'] = read_array(glialDoc.getElementsByTagName("tilt")[0].childNodes[0].nodeValue);
	data['tilt_bin'] = read_array(glialDoc.getElementsByTagName("tilt_bin")[0].childNodes[0].nodeValue);
	data['rotation'] = read_array(glialDoc.getElementsByTagName("rotation")[0].childNodes[0].nodeValue);
	data['rotation_bin'] = read_array(glialDoc.getElementsByTagName("rotation_bin")[0].childNodes[0].nodeValue);
	data['pvt'] = read_2d_array(glialDoc.getElementsByTagName("pvt")[0].childNodes[0].nodeValue);
	data['pvt_pbin'] = read_array(glialDoc.getElementsByTagName("pvt_xbin")[0].childNodes[0].nodeValue);
	data['pvt_tbin'] = read_array(glialDoc.getElementsByTagName("pvt_ybin")[0].childNodes[0].nodeValue);
	data['pvr'] = read_2d_array(glialDoc.getElementsByTagName("pvr")[0].childNodes[0].nodeValue);
	data['pvr_pbin'] = read_array(glialDoc.getElementsByTagName("pvr_xbin")[0].childNodes[0].nodeValue);
	data['pvr_rbin'] = read_array(glialDoc.getElementsByTagName("pvr_ybin")[0].childNodes[0].nodeValue);
	data['tvr'] = read_2d_array(glialDoc.getElementsByTagName("tvr")[0].childNodes[0].nodeValue);
	data['tvr_tbin'] = read_array(glialDoc.getElementsByTagName("tvr_xbin")[0].childNodes[0].nodeValue);
	data['tvr_rbin'] = read_array(glialDoc.getElementsByTagName("tvr_ybin")[0].childNodes[0].nodeValue);
	
	data['pvtvr'] = read_3d_array(glialDoc.getElementsByTagName("pvtvr")[0].textContent);
	data['pvtvr_pbin'] = read_array(glialDoc.getElementsByTagName("pvtvr_xbin")[0].childNodes[0].nodeValue);
	data['pvtvr_tbin'] = read_array(glialDoc.getElementsByTagName("pvtvr_ybin")[0].childNodes[0].nodeValue);
	data['pvtvr_rbin'] = read_array(glialDoc.getElementsByTagName("pvtvr_zbin")[0].childNodes[0].nodeValue);
	//alert(data['position']);
	
	if (document.getElementById("logbinSelect").value == "true") {
		logbin = true;
		} else {
		logbin = false;
		}
	
	if (document.getElementById("sidechainSelect").value == "true") {
		sidechain = true;
		} else {
		sidechain = false;
		}
	
	if (document.getElementById("rotateViewSelect").value == "true") {
		rotateView = true;
		} else {
		rotateView = false;
		}
	
	drawHistogram(data['position'],data['position_bin'], 'position',-1,false,document.getElementById("colorbySelect").value,logbin)
	drawHistogram(data['tilt'],data['tilt_bin'], 'tilt',-1,false,document.getElementById("colorbySelect").value,logbin)
	drawHistogram(data['rotation'],data['rotation_bin'], 'rotation',-1,false,document.getElementById("colorbySelect").value,logbin)
       //alert(document.getElementById("datasetSelect").value)
    draw2dHistogram(data['pvt'],data['pvt_pbin'],data['pvt_tbin'],"pvt",-1,false,false,logbin);
    draw2dHistogram(data['pvr'],data['pvr_pbin'],data['pvr_rbin'],"pvr",-1,false,false,logbin);
    draw2dHistogram(data['tvr'],data['tvr_tbin'],data['tvr_rbin'],"tvr",-1,false,false,logbin);
    tm_modal_position = get_modal_value(data['position'],data['position_bin']);
    tm_modal_tilt = get_modal_value(data['tilt'],data['tilt_bin']);
    tm_modal_rotation = get_modal_value(data['rotation'],data['rotation_bin']);
    
    //modal_orientation = get_3d_modal_value (data['pvtvr'], data['pvtvr_pbin'],data['pvtvr_tbin'],data['pvtvr_rbin'])
    //Nd modal orientation is a bit deceptive; look for position modal orientation, then find its modal tilt, then rotation
    modal_orientation = get_modal_value (data['pvt'][tm_modal_position[1]], data['pvt_tbin']);
    modal_rotation = get_modal_value(data['pvtvr'][tm_modal_position[1]][modal_orientation[1]],data['pvtvr_rbin']);
    
	//helix_cartoon("transmembrane_cartoon",tm_modal_position,tm_modal_tilt,360-tm_modal_rotation,sidechain,data['sequence'],rotateView);
	helix_cartoon("transmembrane_cartoon",tm_modal_position[0],modal_orientation[0],360-modal_rotation[0],sidechain,data['sequence'],rotateView,data['thickness']);
	document.getElementById("header").innerHTML= "<h1>Sidekick Detail-" + data['sequence'] + "</h1>"
}

function draw_cartoons(array_of_identities){
	for (id in array_of_identities) {
		try {
		goofy_cartoon(array_of_identities[id]);
		}
		catch (e) {
		}
	}
}

function goofy_cartoon(whichDiv){
	//document.getElementById('canvas').width=400;
	//document.getElementById('canvas').height=300;
	var data = new Array();
	
	xmlglial=new XMLHttpRequest();
	xmlglial.open("GET",whichDiv+"/"+document.getElementById("datasetSelect").value+"?",false);
	xmlglial.send();
	glialDoc=xmlglial.responseXML;
	data['sequence'] = glialDoc.getElementsByTagName("sequence")[0].childNodes[0].nodeValue
	data['thickness'] = parseInt(glialDoc.getElementsByTagName("thickness")[0].childNodes[0].nodeValue)
	
	data['position'] = read_array(glialDoc.getElementsByTagName("position")[0].childNodes[0].nodeValue);
	data['position_bin'] = read_array(glialDoc.getElementsByTagName("position_bin")[0].childNodes[0].nodeValue);
	data['tilt'] = read_array(glialDoc.getElementsByTagName("tilt")[0].childNodes[0].nodeValue);
	data['tilt_bin'] = read_array(glialDoc.getElementsByTagName("tilt_bin")[0].childNodes[0].nodeValue);
	data['rotation'] = read_array(glialDoc.getElementsByTagName("rotation")[0].childNodes[0].nodeValue);
	data['rotation_bin'] = read_array(glialDoc.getElementsByTagName("rotation_bin")[0].childNodes[0].nodeValue);
	data['pvt'] = read_2d_array(glialDoc.getElementsByTagName("pvt")[0].childNodes[0].nodeValue);
	data['pvt_pbin'] = read_array(glialDoc.getElementsByTagName("pvt_xbin")[0].childNodes[0].nodeValue);
	data['pvt_tbin'] = read_array(glialDoc.getElementsByTagName("pvt_ybin")[0].childNodes[0].nodeValue);
	data['pvr'] = read_2d_array(glialDoc.getElementsByTagName("pvr")[0].childNodes[0].nodeValue);
	data['pvr_pbin'] = read_array(glialDoc.getElementsByTagName("pvr_xbin")[0].childNodes[0].nodeValue);
	data['pvr_rbin'] = read_array(glialDoc.getElementsByTagName("pvr_ybin")[0].childNodes[0].nodeValue);
	data['tvr'] = read_2d_array(glialDoc.getElementsByTagName("tvr")[0].childNodes[0].nodeValue);
	data['tvr_tbin'] = read_array(glialDoc.getElementsByTagName("tvr_xbin")[0].childNodes[0].nodeValue);
	data['tvr_rbin'] = read_array(glialDoc.getElementsByTagName("tvr_ybin")[0].childNodes[0].nodeValue);
	data['pvtvr'] = read_3d_array(glialDoc.getElementsByTagName("pvtvr")[0].textContent);
	data['pvtvr_pbin'] = read_array(glialDoc.getElementsByTagName("pvtvr_xbin")[0].childNodes[0].nodeValue);
	data['pvtvr_tbin'] = read_array(glialDoc.getElementsByTagName("pvtvr_ybin")[0].childNodes[0].nodeValue);
	data['pvtvr_rbin'] = read_array(glialDoc.getElementsByTagName("pvtvr_zbin")[0].childNodes[0].nodeValue);
	//alert(data['position']);
	if (document.getElementById("rotateViewSelect").value == "true") {
		rotateView = true;
		} else {
		rotateView = false;
		}
    tm_modal_position = get_modal_value(data['position'],data['position_bin']);
    modal_orientation = get_modal_value (data['pvt'][tm_modal_position[1]], data['pvt_tbin']);
    modal_rotation = get_modal_value(data['pvtvr'][tm_modal_position[1]][modal_orientation[1]],data['pvtvr_rbin']);
	helix_cartoon(whichDiv,tm_modal_position[0],modal_orientation[0],360-modal_rotation[0],true,data['sequence'],rotateView,data['thickness']);
}
