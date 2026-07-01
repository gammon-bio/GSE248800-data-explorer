// idle-timeout.js — after IDLE_MINUTES of no interaction, close the Shiny
// WebSocket so a forgotten open tab stops holding the server awake. Once no
// session holds a connection, a serverless host (Railway) can sleep the
// container. The user reloads the page to resume (cold start reconnects).
(function () {
  var mins = window.IDLE_MINUTES || 15;
  var IDLE_MS = mins * 60 * 1000;
  var timer = null;

  function pause() {
    try {
      if (window.Shiny && Shiny.shinyapp && Shiny.shinyapp.$socket) {
        Shiny.shinyapp.$socket.close();   // no auto-reconnect (allowReconnect not set)
      }
    } catch (e) { /* ignore */ }
    var o = document.getElementById("idle-overlay");
    if (o) o.style.display = "flex";
  }

  function reset() {
    if (timer) clearTimeout(timer);
    timer = setTimeout(pause, IDLE_MS);
  }

  ["mousemove", "mousedown", "keydown", "scroll", "touchstart", "click", "wheel"]
    .forEach(function (ev) { document.addEventListener(ev, reset, { passive: true }); });

  if (document.readyState !== "loading") reset();
  else document.addEventListener("DOMContentLoaded", reset);
})();
