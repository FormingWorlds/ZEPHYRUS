/* Render mermaid diagrams written as fenced code blocks.
 *
 * The site builder emits a ```mermaid fence as a highlighted code block, so
 * the diagram source reaches the page as text. This finds those blocks,
 * replaces each with the element mermaid renders into, and renders it in the
 * colors of the active scheme. Toggling light and dark mode renders every
 * diagram again from its stored source. Without the library (offline, or a
 * blocked CDN) the page keeps the code block and stays readable.
 *
 * Diagrams assign node classes and leave the colors to this file: `regime`
 * for an outcome, `decision` for a test, `stage` for an intermediate step.
 */
(function () {
  var KEYWORDS = /^\s*(flowchart|graph|sequenceDiagram|classDiagram|stateDiagram(-v2)?|erDiagram|journey|gantt|pie|mindmap|timeline|quadrantChart|gitGraph)\b/;

  // One palette per scheme, on the PROTEUS colors of stylesheets/extra.css.
  var PALETTES = {
    light: {
      vars: {
        background: '#F2F5F7',
        mainBkg: '#FFFFFF',
        primaryColor: '#FFFFFF',
        primaryTextColor: '#0E131B',
        primaryBorderColor: '#8A97A6',
        lineColor: '#4A5868',
        textColor: '#0E131B',
        edgeLabelBackground: '#F2F5F7',
      },
      classes: {
        regime: 'fill:#1B6FA8,stroke:#14406B,color:#FFFFFF',
        decision: 'fill:#FFFFFF,stroke:#8A97A6,color:#0E131B',
        stage: 'fill:#E3EEF6,stroke:#1B6FA8,color:#0E131B',
      },
    },
    dark: {
      vars: {
        background: '#05070B',
        mainBkg: '#0E131B',
        primaryColor: '#0E131B',
        primaryTextColor: '#E6EDF3',
        primaryBorderColor: '#4FA3D9',
        lineColor: '#A8B4C2',
        textColor: '#E6EDF3',
        edgeLabelBackground: '#05070B',
      },
      classes: {
        regime: 'fill:#1B6FA8,stroke:#A8D4E8,color:#FFFFFF',
        decision: 'fill:#0E131B,stroke:#8A97A6,color:#E6EDF3',
        stage: 'fill:#0E2A45,stroke:#4FA3D9,color:#E6EDF3',
      },
    },
  };

  var targets = [];

  function collect() {
    // Both selectors can match inside one block, so key on the element
    // that gets replaced and keep the first source seen for it.
    var seen = new Map();
    document.querySelectorAll('pre code, div.highlight pre').forEach(function (node) {
      var text = node.textContent || '';
      if (!KEYWORDS.test(text)) {
        return;
      }
      var host = node.closest('div.highlight') || node.closest('pre');
      if (host && !seen.has(host)) {
        seen.set(host, text);
      }
    });
    return Array.from(seen, function (entry) {
      return { host: entry[0], source: entry[1] };
    });
  }

  function darkMode() {
    var scheme = document.body.getAttribute('data-md-color-scheme');
    if (scheme) {
      return scheme !== 'default';
    }
    return window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches;
  }

  function themed(source, palette) {
    // Class colors come from the palette, so a diagram's own classDef lines
    // for these names are replaced rather than combined.
    var names = Object.keys(palette.classes);
    var lines = source.split('\n').filter(function (line) {
      var match = line.match(/^\s*classDef\s+(\w+)/);
      return !(match && names.indexOf(match[1]) >= 0);
    });
    names.forEach(function (name) {
      lines.push('    classDef ' + name + ' ' + palette.classes[name]);
    });
    return lines.join('\n');
  }

  function render() {
    var palette = darkMode() ? PALETTES.dark : PALETTES.light;
    window.mermaid.initialize({
      startOnLoad: false,
      securityLevel: 'strict',
      theme: 'base',
      themeVariables: palette.vars,
      flowchart: { htmlLabels: true, useMaxWidth: true },
    });
    targets.forEach(function (target) {
      target.removeAttribute('data-processed');
      target.textContent = themed(target.dataset.source, palette);
    });
    window.mermaid.run({ nodes: targets });
  }

  function setup() {
    if (typeof window.mermaid === 'undefined') {
      return;
    }
    targets = collect().map(function (block) {
      var target = document.createElement('div');
      target.className = 'mermaid';
      target.dataset.source = block.source;
      block.host.replaceWith(target);
      return target;
    });
    if (!targets.length) {
      return;
    }
    render();
    new MutationObserver(render).observe(document.body, {
      attributes: true,
      attributeFilter: ['data-md-color-scheme'],
    });
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', setup);
  } else {
    setup();
  }
})();
