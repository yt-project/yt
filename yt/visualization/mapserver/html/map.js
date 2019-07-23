function setFullScreen () {
    $("#map").width($(window).width());
    $("#map").height($(window).height());
}

var SearchWidget = function () {
  var obj = {
    filter: function (searchStrs) {
      console.log("filtering on " + searchStrs);
      this._selector.each(function(i, el) {
        var val = $(el).text();
        // Search
        var matched = searchStrs.map((str) => {
          return val.indexOf(str) !== -1;
        }).reduce((reduced, result) => {
          return reduced && result;
        }, true);
        if (matched) {
          $(el).show();
        } else {
          $(el).hide();
        }
      });
    },
    init: function () {
      var self = this;

      var searchElement = $('<div id="filter"><input type="text" placeholder="Filter layers"></div>');
      var selector = $('.leaflet-control-layers-list label');

      this._selector = selector;

      // Add input in the DOM
      selector.first().parent().prepend(searchElement);

      // Listen to keyboard input
      $('#filter input').keyup(function(ev) {
        const val = $(this).val();
        self.filter(val.split(" "));
      });

    },
    _selector: null
  };
  obj.init();
  return obj;
};
$(document).ready(function() {
  // Initialize to full screen
  setFullScreen();
  // initialize the map on the "map" div with a given center and zoom
  $.getJSON('/list', function(data) {
    var layers = [],
        layer_groups = [],
        default_layer = [null];
    var layer_group = {};

    // Loop over field types
    for (var type in data['data']) {
      var dtype = data['data'][type];

      // Loop over fields of given type
      for (var field in dtype) {
        var loc = dtype[field]
        var field = loc[0],
            active = loc[1],
            url = 'map/' + field[0] + ',' + field[1] + '/{z}/{x}/{y}.png';

        // Create new layer
        var layer = new L.TileLayer(url, {id: 'MapID', maxzoom: 18});

        // Create readable name
        human_name = field.join(' ');

        // Store it
        layers.push(layer);
        layer_group[human_name] = layer;
        if (active) {
          default_layer[0] = layer;
        }
      }
    }
    var map = new L.Map('map', {
      crs: L.CRS.Simple,
      center: new L.LatLng(-128, -128),
      zoom: 4,
      layers: default_layer
    });

    L.control.layers(layer_group).addTo(map);

    var unit = data['unit'], px2unit = data['px2unit'], decimals = 2;
    var fmt = (n) => {
      return L.NumberFormatter.round(n, decimals, ".")
    };
    L.control.coordinates({
      position: "bottomleft", //optional default "bootomright"
      decimals: 2, //optional default 4
      decimalSeperator: ".", //optional default "."
      enableUserInput: false, //optional default true
      useDMS: false, //optional default false
      useLatLngOrder: false, //ordering of labels, default false-> lng-lat
      markerType: L.marker, //optional default L.marker
      labelFormatterLng : (lng) => {
        return fmt((lng+128)*px2unit) + " " + unit
      }, //optional default none,
      labelFormatterLat : (lat) => {
        return fmt((lat+128)*px2unit) + " " + unit
      }, //optional default none
    }).addTo(map);

    // Search widget
    var search = SearchWidget();
  });

  // Resize map automatically
  $(window).resize(setFullScreen);
});
