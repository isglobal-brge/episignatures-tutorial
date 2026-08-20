(function () {
  const PAGE_KEY = "episig-page-theme";
  const CODE_KEY = "episig-code-theme";
  const COLLAPSE_KEY = "episig-theme-switcher-collapsed";
  const root = document.documentElement;

  function normalize(value) {
    return value === "dark" ? "dark" : "light";
  }

  function setTheme(kind, value) {
    const clean = normalize(value);
    if (kind === "page") {
      root.setAttribute("data-page-theme", clean);
      localStorage.setItem(PAGE_KEY, clean);
      const meta = document.querySelector('meta[name="theme-color"]');
      if (meta) meta.setAttribute("content", clean === "dark" ? "#0f172a" : "#0f766e");
    }
    if (kind === "code") {
      root.setAttribute("data-code-theme", clean);
      localStorage.setItem(CODE_KEY, clean);
    }
    updateButtons();
  }

  function updateButtons() {
    document.querySelectorAll("[data-epi-theme-kind]").forEach(function (button) {
      const kind = button.getAttribute("data-epi-theme-kind");
      const value = button.getAttribute("data-epi-theme-value");
      const current = root.getAttribute("data-" + kind + "-theme") || "light";
      button.setAttribute("aria-pressed", String(current === value));
    });
  }

  function buildSwitcher() {
    const switcher = document.createElement("aside");
    switcher.className = "epi-theme-switcher";
    switcher.setAttribute("aria-label", "Display theme controls");
    switcher.innerHTML = `
      <div class="epi-theme-switcher__title">
        <span>Display</span>
        <button type="button" class="epi-theme-switcher__close" aria-label="Collapse display controls">-</button>
      </div>
      <div class="epi-theme-switcher__body">
        <div class="epi-theme-switcher__row">
          <div class="epi-theme-switcher__label">Text</div>
          <div class="epi-theme-switcher__buttons">
            <button type="button" data-epi-theme-kind="page" data-epi-theme-value="light">Light</button>
            <button type="button" data-epi-theme-kind="page" data-epi-theme-value="dark">Dark</button>
          </div>
        </div>
        <div class="epi-theme-switcher__row">
          <div class="epi-theme-switcher__label">Code</div>
          <div class="epi-theme-switcher__buttons">
            <button type="button" data-epi-theme-kind="code" data-epi-theme-value="light">Light</button>
            <button type="button" data-epi-theme-kind="code" data-epi-theme-value="dark">Dark</button>
          </div>
        </div>
      </div>
    `;

    switcher.addEventListener("click", function (event) {
      const button = event.target.closest("button");
      if (!button) return;
      if (button.classList.contains("epi-theme-switcher__close")) {
        switcher.classList.toggle("is-collapsed");
        const collapsed = switcher.classList.contains("is-collapsed");
        localStorage.setItem(COLLAPSE_KEY, collapsed ? "true" : "false");
        button.textContent = collapsed ? "Aa" : "-";
        button.setAttribute("aria-label", collapsed ? "Open display controls" : "Collapse display controls");
        return;
      }
      const kind = button.getAttribute("data-epi-theme-kind");
      const value = button.getAttribute("data-epi-theme-value");
      if (kind && value) setTheme(kind, value);
    });

    if (localStorage.getItem(COLLAPSE_KEY) === "true") {
      switcher.classList.add("is-collapsed");
      const close = switcher.querySelector(".epi-theme-switcher__close");
      close.textContent = "Aa";
      close.setAttribute("aria-label", "Open display controls");
    }

    document.body.appendChild(switcher);
    updateButtons();
  }

  setTheme("page", localStorage.getItem(PAGE_KEY) || "light");
  setTheme("code", localStorage.getItem(CODE_KEY) || "light");

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", buildSwitcher);
  } else {
    buildSwitcher();
  }
})();
