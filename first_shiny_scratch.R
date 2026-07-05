library(shiny)
library(htmltools)

ui <- fluidPage(
  titlePanel("Dynamic MapGL 3D Globe to 2D Leaflet Transition"),
  tags$head(
    # Upgraded to MapLibre GL v5.8.0 for native 3D Globe Support
    tags$link(rel = "stylesheet", href = "https://unpkg.com/maplibre-gl@5.8.0/dist/maplibre-gl.css"),
    tags$script(src = "https://unpkg.com/maplibre-gl@5.8.0/dist/maplibre-gl.js"),
    
    # Load Leaflet CSS & JS
    tags$link(rel = "stylesheet", href = "https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"),
    tags$script(src = "https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"),
    
    # Styling to handle smooth overlay swap
    tags$style(HTML("
      #map-container {
        position: relative;
        width: 100%;
        height: 700px;
        background: #0d1117;
        border-radius: 8px;
        overflow: hidden;
      }
      .map-layer {
        position: absolute;
        top: 0;
        left: 0;
        width: 100%;
        height: 100%;
        transition: opacity 0.3s ease-in-out;
      }
      .hidden {
        opacity: 0 !important;
        pointer-events: none;
      }
    "))
  ),
  
  sidebarLayout(
    sidebarPanel(
      p("Zoom into the globe to seamlessly transition into a 2D Leaflet street map."),
      textOutput("zoom_ui"),
      br(),
      textOutput("click_ui")
    ),
    mainPanel(
      div(id = "map-container",
          div(id = "mapgl", class = "map-layer"),
          div(id = "leaflet", class = "map-layer hidden")
      ),
      
      tags$script(HTML("
        document.addEventListener('DOMContentLoaded', function() {
          const ZOOM_THRESHOLD = 4; 
          let currentMode = 'globe'; 
          
          // 1. Initialize MapLibre GL Map with a Terrain/Satellite Style
          const mapgl = new maplibregl.Map({
            container: 'mapgl',
            style: {
              'version': 8,
              'sources': {
                'satellite-terrain': {
                  'type': 'raster',
                  'tiles': ['https://tiles.maps.eox.at/wmts/1.0.0/s2cloudless-2020_3857/default/g/{z}/{y}/{x}.jpg'],
                  'tileSize': 256,
                  'attribution': '© EOX IT Services GmbH'
                }
              },
              'layers': [
                {
                  'id': 'satellite-layer',
                  'type': 'raster',
                  'source': 'satellite-terrain'
                }
              ]
            },
            center: [18.56, 54.44], 
            zoom: 1.5,
            maxPitch: 85
          });

          // 2. Explicitly apply Globe projection once the basemap style loads
          mapgl.on('style.load', () => {
            mapgl.setProjection({ type: 'globe' });
          });

          // 3. Initialize Leaflet instance (hidden initially)
          const leaflet = L.map('leaflet').setView([54.44, 18.56], 5);
          L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
            attribution: '© OpenStreetMap contributors'
          }).addTo(leaflet);

          // Synchronize MapGL coordinates -> Leaflet
          function syncToLeaflet() {
            const center = mapgl.getCenter();
            const zoom = mapgl.getZoom();
            leaflet.setView([center.lat, center.lng], Math.round(zoom) + 1, { animate: false });
          }

          // Synchronize Leaflet coordinates -> MapGL
          function syncToMapGL() {
            const center = leaflet.getCenter();
            const zoom = leaflet.getZoom();
            mapgl.setCenter([center.lng, center.lat]);
            mapgl.setZoom(zoom - 1);
          }

          // Capture Location Selection on MapGL (Globe)
          mapgl.on('click', (e) => {
            if (currentMode === 'globe') {
              Shiny.setInputValue('selected_location', {
                lat: e.lngLat.lat,
                lng: e.lngLat.lng,
                engine: 'MapLibre Globe'
              }, {priority: 'event'});
            }
          });

          // Capture Location Selection on Leaflet (2D Map)
          leaflet.on('click', (e) => {
            if (currentMode === 'leaflet') {
              Shiny.setInputValue('selected_location', {
                lat: e.latlng.lat,
                lng: e.latlng.lng,
                engine: 'Leaflet 2D'
              }, {priority: 'event'});
            }
          });

          // Monitor MapLibre GL Zooms
          mapgl.on('zoom', () => {
            const zoom = mapgl.getZoom();
            Shiny.setInputValue('current_zoom', zoom.toFixed(2));
            
            if (zoom >= ZOOM_THRESHOLD && currentMode === 'globe') {
              currentMode = 'leaflet';
              syncToLeaflet();
              document.getElementById('leaflet').classList.remove('hidden');
              document.getElementById('mapgl').classList.add('hidden');
              setTimeout(() => leaflet.invalidateSize(), 50); 
            }
          });

          // Monitor Leaflet zooms/pans to drop back to globe when zooming out
          leaflet.on('zoomend moveend', () => {
            if (currentMode === 'leaflet') {
              const zoom = leaflet.getZoom();
              Shiny.setInputValue('current_zoom', (zoom - 1).toFixed(2));
              
              if (zoom < (ZOOM_THRESHOLD + 1)) {
                currentMode = 'globe';
                syncToMapGL();
                document.getElementById('mapgl').classList.remove('hidden');
                document.getElementById('leaflet').classList.add('hidden');
              }
            }
          });
        });
      "))
    )
  )
)

server <- function(input, output, session) {
  output$zoom_ui <- renderText({
    req(input$current_zoom)
    paste("Current MapGL Zoom Level:", input$current_zoom)
  })
  
  output$click_ui <- renderText({
    req(input$selected_location)
    loc <- input$selected_location
    paste0("Selected via ", loc$engine, ": [Lat: ", round(loc$lat, 4), ", Lng: ", round(loc$lng, 4), "]")
  })
}

shinyApp(ui, server)