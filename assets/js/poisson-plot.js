/*
 * Interactive render of the DOLFINx Poisson demo solution.
 *
 * The mesh (561 vertices, 1024 triangles) lives in assets/plots/poisson.json.
 * Plotly's gl3d bundle is 1.6MB, so neither it nor the mesh is fetched until
 * the plot is close to the viewport.
 */
(function () {
  "use strict";

  var PLOTLY_SRC = "https://cdn.plot.ly/plotly-gl3d-2.20.0.min.js";
  var PLOTLY_SRI = "sha384-8ordnAyND+ShSzrlsx3CJMv3wFJgi7oh8eNgyKmc9Lu6qV2mY9QRBT4ROJAr351E";
  var MESH_SRC = "/assets/plots/poisson.json";

  /* Closer than Plotly's default eye, which is set up to frame a unit cube;
   * this domain is a flat 2 x 1 sheet and would otherwise sit small in frame. */
  var CAM = { x: 1.22, y: -1.58, z: 1.02 };

  var holder = document.getElementById("poisson-plot");
  if (!holder) return;

  var status = document.getElementById("poisson-status");
  var controls = document.getElementById("poisson-controls");

  function fail(msg) {
    if (status) status.textContent = msg;
  }

  function loadPlotly() {
    return new Promise(function (resolve, reject) {
      if (window.Plotly) return resolve();
      var s = document.createElement("script");
      s.src = PLOTLY_SRC;
      s.integrity = PLOTLY_SRI;
      s.crossOrigin = "anonymous";
      s.onload = resolve;
      s.onerror = function () { reject(new Error("Plotly failed to load")); };
      document.head.appendChild(s);
    });
  }

  /* Unique edges of the triangulation, as a single scatter3d polyline broken
   * by nulls. Derived here rather than stored, since it costs nothing at this
   * mesh size and would double the payload on disk.
   *
   * The edges are lifted by a fraction of the z range: drawn exactly on the
   * surface they z-fight with it and vanish in the "Both" view. */
  function wireframe(m) {
    var lift = 0.012 * (Math.max.apply(null, m.z) - Math.min.apply(null, m.z));
    var seen = Object.create(null);
    var ex = [], ey = [], ez = [];
    function edge(a, b) {
      var key = a < b ? a + "_" + b : b + "_" + a;
      if (seen[key]) return;
      seen[key] = true;
      ex.push(m.x[a], m.x[b], null);
      ey.push(m.y[a], m.y[b], null);
      ez.push(m.z[a] + lift, m.z[b] + lift, null);
    }
    for (var t = 0; t < m.i.length; t++) {
      edge(m.i[t], m.j[t]);
      edge(m.j[t], m.k[t]);
      edge(m.k[t], m.i[t]);
    }
    return { x: ex, y: ey, z: ez };
  }

  function render(m) {
    var w = wireframe(m);

    var surface = {
      type: "mesh3d",
      x: m.x, y: m.y, z: m.z,
      i: m.i, j: m.j, k: m.k,
      intensity: m.z,
      colorscale: "Viridis",
      flatshading: false,
      /* Lighting is kept deliberately flat. Plotly derives mesh3d normals per
       * facet, so a strong diffuse term makes the triangulation read as a
       * checkerboard; letting the colourmap carry the shape avoids that. */
      lighting: { ambient: 0.82, diffuse: 0.32, specular: 0.04, roughness: 0.9, fresnel: 0.05 },
      lightposition: { x: 1.4, y: -1.6, z: 2.4 },
      colorbar: {
        title: { text: "u<sub>h</sub>", side: "right" },
        thickness: 14,
        len: 0.7,
        outlinewidth: 0,
        tickfont: { size: 11 }
      },
      hovertemplate: "x = %{x:.3f}<br>y = %{y:.3f}<br>u<sub>h</sub> = %{z:.4f}<extra></extra>",
      name: "solution"
    };

    var mesh = {
      type: "scatter3d",
      mode: "lines",
      x: w.x, y: w.y, z: w.z,
      line: { color: "rgba(30,30,30,0.55)", width: 1 },
      hoverinfo: "skip",
      visible: false,
      name: "mesh"
    };

    var layout = {
      /* The right margin reserves the colourbar's space so the scene does not
       * jump sideways when the surface is toggled off. */
      margin: { l: 0, r: 70, t: 0, b: 0 },
      paper_bgcolor: "rgba(0,0,0,0)",
      showlegend: false,
      scene: {
        aspectmode: "data",
        camera: { eye: startEye() },
        xaxis: { title: "x", gridcolor: "#ddd", zerolinecolor: "#bbb", showbackground: false },
        yaxis: { title: "y", gridcolor: "#ddd", zerolinecolor: "#bbb", showbackground: false },
        zaxis: { title: "u<sub>h</sub>", gridcolor: "#ddd", zerolinecolor: "#bbb", showbackground: false }
      }
    };

    var config = {
      responsive: true,
      displaylogo: false,
      modeBarButtonsToRemove: ["toImage", "resetCameraLastSave3d"]
    };

    return window.Plotly.newPlot(holder, [surface, mesh], layout, config).then(function () {
      if (status) status.remove();
      if (controls) controls.hidden = false;
      wireControls();
      spin();
    });
  }

  function wireControls() {
    if (!controls) return;
    /* [surface visible, wireframe visible] per view */
    var views = { surface: [true, false], mesh: [false, true], both: [true, true] };
    controls.addEventListener("click", function (ev) {
      var btn = ev.target.closest("button[data-view]");
      if (!btn) return;
      var v = views[btn.dataset.view];
      if (!v) return;
      window.Plotly.restyle(holder, { visible: v }, [0, 1]);
      Array.prototype.forEach.call(controls.querySelectorAll("button[data-view]"), function (b) {
        b.setAttribute("aria-pressed", String(b === btn));
      });
    });
  }

  var SWING = 0.32;                                   /* radians of opening orbit */
  var TARGET = Math.atan2(CAM.y, CAM.x);
  var RADIUS = Math.sqrt(CAM.x * CAM.x + CAM.y * CAM.y);
  var reduceMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  function eyeAt(angle) {
    return { x: RADIUS * Math.cos(angle), y: RADIUS * Math.sin(angle), z: CAM.z };
  }

  /* The opening orbit runs *into* CAM rather than away from it, so the plot
   * comes to rest on the framing that was tuned, not 20 degrees past it. */
  function startEye() {
    return reduceMotion ? eyeAt(TARGET) : eyeAt(TARGET - SWING);
  }

  /* A slow orbit on arrival, so it reads as interactive. Any input stops it
   * for good — nothing is more irritating than a plot that fights the cursor. */
  function spin() {
    if (reduceMotion) return;

    var angle = TARGET - SWING;
    var running = true;
    var timer = null;

    function stop() {
      running = false;
      if (timer) clearInterval(timer);
      holder.removeEventListener("pointerdown", stop);
      holder.removeEventListener("wheel", stop);
    }
    holder.addEventListener("pointerdown", stop, { once: true });
    holder.addEventListener("wheel", stop, { once: true, passive: true });

    timer = setInterval(function () {
      if (!running) return;
      angle = Math.min(angle + 0.0022, TARGET);
      window.Plotly.relayout(holder, { "scene.camera.eye": eyeAt(angle) });
      if (angle >= TARGET) stop();
    }, 40);
  }

  function boot() {
    fail("Loading…");
    loadPlotly()
      .then(function () { return fetch(MESH_SRC); })
      .then(function (r) {
        if (!r.ok) throw new Error("mesh " + r.status);
        return r.json();
      })
      .then(render)
      .catch(function (e) {
        fail("The interactive plot could not be loaded (" + e.message + ").");
      });
  }

  if (!("IntersectionObserver" in window)) {
    boot();
    return;
  }

  var io = new IntersectionObserver(function (entries) {
    if (entries.some(function (e) { return e.isIntersecting; })) {
      io.disconnect();
      boot();
    }
  }, { rootMargin: "300px" });
  io.observe(holder);
})();
