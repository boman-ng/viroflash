#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntegrationStatus {
    NotObserved,
    DiagnosticEvidenceObserved,
}

impl IntegrationStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::NotObserved => "NOT_OBSERVED",
            Self::DiagnosticEvidenceObserved => "DIAGNOSTIC_EVIDENCE_OBSERVED",
        }
    }
}
