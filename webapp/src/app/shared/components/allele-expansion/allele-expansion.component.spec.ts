import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { AlleleExpansionComponent } from './allele-expansion.component';

describe('AlleleExpansionComponent', () => {
  let component: AlleleExpansionComponent;
  let fixture: ComponentFixture<AlleleExpansionComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ AlleleExpansionComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(AlleleExpansionComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
