function wire() {
  const homepage = "https://proteus-framework.org/";

  const logo = document.querySelector(".md-header__button.md-logo");
  if (logo) logo.href = homepage;

  const title = document.querySelector(".md-header__title[data-md-component='header-title']");
  if (title && !title.dataset.spiderWired) {
    title.dataset.spiderWired = "1";
    title.style.cursor = "pointer";

    // derive the hosted docs base path dynamically, else "/" (mkdocs serve)
    const pathMatch = location.pathname.match(/^\/([A-Z0-9_-]+)(?:\/|$)/);
    const docsHome = pathMatch
      ? location.origin + "/" + pathMatch[1] + "/"
      : location.origin + "/";

    title.addEventListener("click", (e) => {
      if (e.target.closest("a, button, input, label")) return;
      window.location.assign(docsHome);
    }, true);
  }
}

document.addEventListener("DOMContentLoaded", wire);
if (window.document$?.subscribe) window.document$.subscribe(wire);
