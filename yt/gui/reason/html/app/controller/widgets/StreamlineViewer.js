/**********************************************************************
The Streamline Viewer Widget

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/



var WidgetStreamlineViewer = function(python_varname, widget_data) {
    this.id = python_varname;
    this.widget_data = widget_data;
    examine = "canvas_" + python_varname;
    var StreamlineViewerStart = function() {
        this.curX = 0;
        this.curY = 0;
        this.dist = 0;
        function updateBasedOnOffset(camera, offset){
        camera.position.x = camera.target.x + offset.x;
        camera.position.y = camera.target.y + offset.y;
        camera.position.z = camera.target.z + offset.z;
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
            program: {
            from: 'ids',
                vs: 'sl-shader-vs',
                fs: 'sl-shader-fs'
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
                c.update();
                }
    
            },
                onMouseWheel: function(e){
                e.stop();
                var offset = PhiloGL.Vec3.scale(camGetOffset(this.camera),
                                1.0 - e.wheel/10.);
                updateBasedOnOffset(this.camera, offset);
                this.camera.update();
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

	    gl.viewport(0, 0, canvas.width, canvas.height);
	    gl.clearColor(0, 0, 0, 1);
	    gl.clearDepth(1);
	    gl.blendFunc(gl.SRC_ALPHA, gl.ONE);
	    gl.enable(gl.BLEND);
	    gl.disable(gl.DEPTH_TEST);
	    program.setUniform('alpha', 0.8);

	    program.setBuffers({

		    'shapeset': {
			attribute: 'aVertexPosition',
			    value: new Float32Array(widget_data['stream_positions']),
			    size: 3
			    },
			
			'shapesetColors': {
			    attribute: 'aVertexColor',
				value: new Float32Array(widget_data['stream_colors']),
				size: 4
			    },
			    });


	    camera.view.id();
	    setInterval(draw, 30/60);
	    var stream_counter =0;
	    //Draw the scene
	    function draw() {
		stream_counter = 0;
		gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
		//Draw Triangle
		program.setUniform('uMVMatrix', camera.view);
		program.setUniform('uPMatrix', camera.projection);
		program.setBuffer('shapeset');
		program.setBuffer('shapesetColors');
		for (var i=0; i<widget_data['n_streamlines']; i++){
		    gl.drawArrays(gl.LINES, stream_counter, widget_data['stream_lengths'][i]-1);
		    gl.drawArrays(gl.LINES, stream_counter+1, widget_data['stream_lengths'][i]-1);
		    stream_counter += widget_data['stream_lengths'][i];
		}
		
	    }
	    }
        });  
    }        

    viewport.get("center-panel").add(
        {
            xtype: 'panel',
            id: "sl_" + python_varname,
            title: "WebGL Streamline Viewer",
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
            listeners: { afterlayout: StreamlineViewerStart },
        }
    );

    viewport.get("center-panel").activate("sl_" + this.id);
    viewport.doLayout();
    this.panel = viewport.get("center-panel").get("sl_" + python_varname);
    this.panel.doLayout();

    this.accept_results = function(payload) { }
}

widget_types['streamline_viewer'] = WidgetStreamlineViewer;
