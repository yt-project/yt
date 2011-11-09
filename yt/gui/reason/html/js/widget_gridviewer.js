/**********************************************************************
The Plot Window Widget

Author: Samuel Skillman <samskillman@gmail.com>
Affiliation: University of Colorado at Boulder
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/

// shim layer with setTimeout fallback
window.requestAnimFrame = (function(){
    return  window.requestAnimationFrame       || 
        window.webkitRequestAnimationFrame || 
        window.mozRequestAnimationFrame    || 
        window.oRequestAnimationFrame      || 
        window.msRequestAnimationFrame     || 
        function(/* function */ callback, /* DOMElement */ element){
            window.setTimeout(callback, 1000 / 60);
        };
})();

var WidgetGridViewer = function(python_varname, widget_data) {
    this.id = python_varname;
    this.widget_data = widget_data;
    examine = "canvas_" + python_varname;
    var draw;
    var GridViewerStart = function() {
        this.curX = 0;
        this.curY = 0;
        this.dist = 0;
        function updateBasedOnOffset(camera, offset){
            camera.position.x = camera.target.x + offset.x;
            camera.position.y = camera.target.y + offset.y;
            camera.position.z = camera.target.z + offset.z;
            draw();
        }
        function camGetOffset(camera){
            return PhiloGL.Vec3.sub(camera.position, camera.target)
        }
        PhiloGL('canvas_' + python_varname, {
            camera: {
		position: {
                    x: 0.5, y: 0.5, z: 5
		},
                target: {
                    x: 0.5, y: 0.5, z: 0.5
                },
            },
            events: {
		onDragStart: function(e) {
                    pos = {
			x: e.x,
			y: e.y
                    };
                    this.curX = e.x;
                    this.curY = e.y;
                    this.dist = camGetOffset(this.camera).norm();
		},
                onDragEnd: function(e) {
                    pos = {
			x: e.x,
			y: e.y
                    };
  		},
                onDragMove: function(e) {
                    var c = this.camera;
                    var off = camGetOffset(c);

                    // Get Horizontal vector
                    var horiz = PhiloGL.Vec3.cross(c.up, 
						   camGetOffset(c))
                    horiz.$scale(1./horiz.norm());

                    if (e.event.button == 0){ // Rotation
			// Do vertical rotation about horizontal vector
			var vert_rot = new PhiloGL.Mat4();
			vert_rot.id();
			vert_rot.$rotateAxis((e.y-this.curY)/100., horiz);
			PhiloGL.Mat4.$mulVec3(vert_rot, off);
			PhiloGL.Mat4.$mulVec3(vert_rot, c.up);
			c.up.$scale(1./c.up.norm());

			// Do horizontal rotation about up vector
			var side_rot = new PhiloGL.Mat4();
			side_rot.id();
			side_rot.$rotateAxis(-(e.x-this.curX)/100., c.up);
			side_rot.$mulVec3(off);
			
			// Update current positions
			this.curX = e.x;
			this.curY = e.y;
			this.dist = off.norm();
			updateBasedOnOffset(c, off);
			this.camera.near = this.dist/100000.0;
			this.camera.far = this.dist*2.0;
			c.update();
                    } else if (e.event.button = 2){ // Right click - transpose
			var tscale = 1.0*off.norm()/512.;
			var move_up = c.up.scale(-(e.y-this.curY)*tscale);
			var move_over = horiz.scale(-(e.x-this.curX)*tscale);
			c.position.$add(move_up);
			c.position.$add(move_over);
			c.target.$add(move_up);
			c.target.$add(move_over);
			// Update current positions
			this.curX = e.x;
			this.curY = e.y;
			this.dist = off.norm();
			this.camera.near = this.dist/100000.0;
			this.camera.far = this.dist*2.0;
			c.update();
                    }
		    draw();
		},
                onMouseWheel: function(e){
                    e.stop();
                    var offset = PhiloGL.Vec3.scale(camGetOffset(this.camera),
						    1.0 - e.wheel/10.);
                    updateBasedOnOffset(this.camera, offset);
		    var dist = offset.norm()
                    this.camera.near = offset.norm()/100000.0;
                    this.camera.far = offset.norm()*2.0;
                    this.camera.update();
                    draw();
		}
            },
            onError: function() {
		alert("An error ocurred while loading the application");
            },
            onLoad: function(app) {
		var gl = app.gl,
                canvas = app.canvas,
                program = app.program,
                scene = app.scene,
                camera = app.camera;
		var grids = new PhiloGL.O3D.Model({
            vertices : widget_data['vertex_positions'],
            drawType : "LINES",
            colors : widget_data['vertex_colors'],
        });
        scene.add(grids);
		gl.viewport(0, 0, canvas.width, canvas.height);
		gl.clearColor(0, 0, 0, 1);

		//examine = camera;
		camera.view.id();
		camera.update();
		
		//Draw the scene
		draw = function() {
	    	    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
                scene.render();
		}

		draw();

	    }
        });  
    }        

    viewport.get("center-panel").add(
        {
            xtype: 'panel',
            id: "gv_" + python_varname,
            title: "WebGL Grid Viewer",
            iconCls: 'graph',
            autoScroll: true,
            layout:'absolute',
            closable: true,
            items: [
                { xtype:'panel',
                  autoEl: {
                      tag: 'canvas',
                      id: 'canvas_' + python_varname,
                      style: 'border: none;',
                      width: 512, height:512
                  },
                  width: 512,
                  height: 512
                }],
            listeners: { afterlayout: GridViewerStart },
        }
    );

    viewport.get("center-panel").activate("gv_" + this.id);
    viewport.doLayout();
    this.panel = viewport.get("center-panel").get("gv_" + python_varname);
    this.panel.doLayout();

    this.accept_results = function(payload) { }
}

widget_types['grid_viewer'] = WidgetGridViewer;
